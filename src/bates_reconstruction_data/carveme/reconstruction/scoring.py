# This file was updated from original version in CarveMe by Adam Bates on 7/7/26

import pandas as pd
import polars as pl


def merge_subunits(genes):
    """ Merge list of protein subunit genes into complex

    Args:
        genes (pandas.Series): list of genes

    Returns:
        str: boolean rule

    """
    genes = genes.dropna()

    if len(genes) == 0:
        return None
    else:
        protein = ' and '.join(sorted(genes))
        if len(genes) > 1:
            return '(' + protein + ')'
        else:
            return protein


def merge_subunits_pl(name):
    genes = pl.col(name).drop_nulls()
    return (
        pl.when(genes.len().gt(1))
        .then("(" + genes.str.join(" and ") + ")")
        .otherwise(genes.first())
    )


def merge_subunit_scores(scores):
    """ Merge scores of all genes in a protein complex.
    Calculates the mean score among all subunits.

    Args:
        scores: individual gene scores

    Returns:
        float: merged score

    """

    return scores.fillna(0).mean()


def merge_subunit_scores_pl(name):
    return pl.col(name).fill_nan(0).mean()


def merge_proteins(proteins):
    """ Merge all isozymes that catalyze a given reaction.
    Automatically removes all isozymes with missing score.

    Args:
        proteins (pandas.Series): list of proteins

    Returns:
        str: boolean rule

    """
    proteins = set(proteins.dropna())
    if not proteins:
        return None

    gpr_str = ' or '.join(sorted(proteins))

    if len(proteins) > 1:
        return '(' + gpr_str + ')'
    else:
        return gpr_str


def merge_proteins_pl(name):
    proteins = pl.col(name).drop_nulls().unique()
    return (
        pl.when(proteins.len().gt(1))
        .then("(" + proteins.str.join(" or ") + ")")
        .otherwise(proteins.first())
    )


def merge_protein_scores(scores):
    """ Merge scores of all isozymes that catalyze a given reaction.
    Calculates the maximum score among all isozymes.

    Args:
        scores (pandas.Series): protein scores

    Returns:
        float: merged score

    """
    return scores.max(skipna=True)


def merge_protein_scores_pl(name):
    return pl.col(name).drop_nans().max()


def reaction_scoring(annotation, gprs, spontaneous_score=0.0, debug_output=None):
    """ Calculate reaction scores using new eggnog output.

    Args:
        annotation (pandas.DataFrame): gene annotation results
        gprs (pandas.DataFrame): BiGG GPR rules
        spontaneous_score (float): score to give to spontaneous reactions (default: 0.0)

    Returns:
        pandas.DataFrame: reaction scores
    """
    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)

    # filter best match for each gene
    gene2gene = annotation.sort_values(by='score', ascending=False) \
                          .groupby('target_gene', as_index=False).apply(lambda x: x.iloc[0])

    # merge with gpr table
    gprs['target_gene'] = gprs.apply(lambda row: f"{row['model']}.{row['gene'][2:]}", axis=1)
    gene_scores = pd.merge(gene2gene, gprs, how='right')

    # add default scores for spontaneous genes
    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT', 'G_PP_s0001'}
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'score'] = spontaneous_score
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'query_gene'] = 'spontaneous'

    # from gene to protein scores
    protein_scores = gene_scores.groupby(['protein', 'reaction', 'model'], as_index=False) \
        .agg({'query_gene': merge_subunits, 'score': merge_subunit_scores})

    protein_scores.rename(columns={'query_gene': 'GPR'}, inplace=True)

    # from protein to reaction scores
    reaction_scores = protein_scores.groupby(['reaction'], as_index=False) \
        .agg({'GPR': merge_proteins, 'score': merge_protein_scores}).dropna()

    avg_score = reaction_scores.query('score > 0')['score'].median()

    if avg_score == 0:
        return None, gene2gene

    reaction_scores['normalized_score'] = (reaction_scores['score'] / avg_score).apply(lambda x: round(x, 1))

    if debug_output:
        gene_scores.to_csv(debug_output + '_gene_scores.tsv', sep='\t', index=False)
        protein_scores.to_csv(debug_output + '_protein_scores.tsv', sep='\t', index=False)
        reaction_scores.to_csv(debug_output + '_reaction_scores.tsv', sep='\t', index=False)

    return reaction_scores, gene2gene

def reaction_scoring_pl(annotation: pd.DataFrame, gprs: pd.DataFrame, spontaneous_score=0.0, debug_output=None):
    gene2gene = (
        pl.from_pandas(annotation)
        .group_by("target_gene")
        .agg(
            pl.col.score.max(),
            pl.col.query_gene.get(pl.col.score.arg_max())
        )
    )

    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT', 'G_PP_s0001'}

    reaction_scores = (
        pl.from_pandas(gprs)
        .with_columns(
            target_gene=pl.col.model + "." + pl.col.gene.str.slice(2)
        )
        .join(gene2gene, on="target_gene", how="left")
        .with_columns(
            score=(
                pl.when(pl.col.gene.is_in(spontaneous))
                .then(spontaneous_score)
                .otherwise(pl.col.score)
            ),
            query_gene=(
                pl.when(pl.col.gene.is_in(spontaneous))
                .then(pl.lit("spontaneous"))
                .otherwise(pl.col.query_gene)
            )
        )
        .group_by("protein", "reaction", "model")
        .agg(
            GPR=merge_subunits_pl("query_gene"),
            score=merge_subunit_scores_pl("score")
        )
        .group_by("reaction")
        .agg(
            GPR=merge_proteins_pl("GPR"),
            score=merge_protein_scores_pl("score")
        )
        .drop_nulls()
    )

    avg_score = (
        reaction_scores
        .filter(pl.col.score.gt(0))
        .get_column("score")
        .median()
    )
    if avg_score is None:
        return None, gene2gene.to_pandas()
    
    reaction_scores = reaction_scores.with_columns(
        normalized_score=(pl.col.score / avg_score).round(1)
    )

    if debug_output:
        reaction_scores.write_csv(debug_output + '_reaction_scores.tsv', separator="\t")

    return reaction_scores.to_pandas(), gene2gene.to_pandas()
