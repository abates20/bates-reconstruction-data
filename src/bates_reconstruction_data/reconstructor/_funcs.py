from typing import Optional, Union, Sequence
from collections import defaultdict
import os
import re

import cobra
from cobra.manipulation.delete import remove_genes
import pandas as pd

from .diamond import Diamond
from . import utils, resources
from ._version import __version__


def run_blast(
    inputfile: Union[str, os.PathLike],
    outputfile: Union[str, os.PathLike],
    database: Union[str, os.PathLike],
    processors: Optional[int] = None,
) -> Union[str, os.PathLike]:
    """
    Runs protein BLAST and saves results.
    """

    if processors is None:
        options = ["--more-sensitive", "--max-target-seqs", "1"]
    else:
        options = ["-p", str(processors), "--more-sensitive", "--max-target-seqs", "1"]

    diamond = Diamond()
    diamond.blastp(database, inputfile, outputfile, *options, capture_output=True)

    return outputfile


def read_blast(blast_output_path: Union[str, os.PathLike]) -> pd.DataFrame:
    """
    Retrieves KEGG hits from blast output.
    """
    kegg_gene_pattern = re.compile(r"([a-z]{3,4}:[\w_]+)")
    default_diamond_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore"
    ]
    hits = (
        pd.read_csv(blast_output_path, sep="\t", names=default_diamond_cols)
        [["qseqid", "sseqid"]]
        .assign(
            sbml_id=lambda x: x.qseqid.apply(utils.sanitize_sbml_id),
            kegg_id=lambda x: x.sseqid.str.extract(kegg_gene_pattern, expand=False)
        )
    )
    return hits


def add_genes_from_organism(blast_hits: pd.DataFrame, organism: str, gene_modelseed: dict[str, list[str]]) -> pd.DataFrame:
    blasted_genes = set(blast_hits.get("kegg_id", []))
    org_genes = _get_org_rxns(gene_modelseed, organism).difference(blasted_genes)
    print(f"Adding {len(org_genes)} from organism {organism}")

    org_df = pd.DataFrame(
        {
            "qseqid": None,
            "sseqid": org_genes,
            "sbml_id": [g.split(":")[-1] for g in org_genes],
            "kegg_id": org_genes
        }
    )
    return pd.concat([blast_hits, org_df], ignore_index=True)


def _get_org_rxns(gene_modelseed: dict[str, list[str]], organism: str) -> set[str]:
    """
    Get genes for organism from reference genome.
    """

    org_genes = []
    for gene in gene_modelseed.keys():
        current = gene.split(":")[0]
        if current == organism:
            org_genes.append(gene)

    return set(org_genes)


def genes_to_rxns(
    gene_df: pd.DataFrame,
    gene_rxn_map: dict[str, list[str]],
) -> defaultdict[str, list[str]]:
    """
    Translates genes to ModelSEED reactions
    """

    rxn_db: defaultdict[str, list[str]] = defaultdict(list)
    for gene_id, sbml_id in gene_df[["sseqid", "sbml_id"]].itertuples(index=False):
        for rxn in gene_rxn_map.get(gene_id, []):
            rxn_db[rxn].append(sbml_id)

    return rxn_db


def create_gene_objects(gene_df: pd.DataFrame, gene_names: dict[str, str]) -> list[cobra.Gene]:
    def _make_gene(row: pd.Series):
        gene = cobra.Gene(row.get("sbml_id"))
        kegg_id = row.get("kegg_id")
        if pd.notna(kegg_id):
            gene.annotation["kegg.genes"] = kegg_id
        if gene.id in gene_names:
            gene.name = gene_names[gene.id].title()
        return gene
    
    genes: pd.Series = gene_df.apply(_make_gene, axis=1)
    return genes.to_list()


def create_model(
    template: Union[cobra.Model, str],
    rxn_db: dict[str, list[str]],
    genes: list[cobra.Gene],
    universal: cobra.Model,
    model_id: Optional[str] = None,
) -> cobra.Model:
    """
    Create draft GENRE and integrate GPRs.
    """
    if isinstance(template, str):
        new_model = resources.get_template_model(template)
        new_model.name = model_id # overwrite the template model ID
    else:
        new_model = template.copy()

    if model_id is not None:
        new_model.id = model_id
    new_model.notes["source"] = f"Reconstructor v{__version__}"

    # Add the genes to the model
    new_model.genes.extend(genes)

    # Add the reactions to the model
    for rxn_id in rxn_db.keys():
        if universal.reactions.has_id(rxn_id):
            rxn: cobra.Reaction = universal.reactions.get_by_id(rxn_id)
            new_model.add_reactions([rxn.copy()])
            gpr = " or ".join(rxn_db[rxn_id])
            new_model.reactions.get_by_id(rxn_id).gene_reaction_rule = gpr

    # Remove any genes not linked to at least one reaction
    remove_genes(new_model, [g for g in new_model.genes if len(g.reactions) == 0])

    # Make sure any extracellular metabolites have exchanges
    utils.add_missing_exchanges(new_model)

    # Set the compartment names using the universal model
    new_model.compartments = universal.compartments

    return new_model


def add_annotation(model: cobra.Model) -> cobra.Model:
    """
    Add gene, metabolite, reaction, and biomass reaction annotations.
    """

    # Genes
    for gene in model.genes:
        gene.annotation["sbo"] = "SBO:0000243"

    # Metabolites
    seed_cpd_pattern = re.compile(r"cpd\d{5}(_\w+)?")
    for cpd in model.metabolites:
        cpd.annotation["sbo"] = "SBO:0000247"
        if seed_cpd_pattern.match(cpd.id):
            cpd.annotation["seed.compound"] = cpd.id.split("_")[0]

    # Reactions
    seed_rxn_pattern = re.compile(r"rxn\d{5}(_\w+)?")
    for rxn in model.reactions:
        if seed_rxn_pattern.match(rxn.id):
            rxn.annotation["seed.reaction"] = rxn.id.split("_")[0]
        compartments = set([x.compartment for x in list(rxn.metabolites)])
        if len(list(rxn.metabolites)) == 1:
            if rxn.id.startswith("EX_"):
                rxn.annotation["sbo"] = "SBO:0000627"  # exchange
            elif rxn.id.startswith("SK_"):
                rxn.annotation["sbo"] = "SBO:0000632" # sink
        elif len(compartments) > 1:
            rxn.annotation["sbo"] = "SBO:0000185"  # transport
        else:
            rxn.annotation["sbo"] = "SBO:0000176"  # metabolic

    return model


def check_model(
    pre_reactions: Sequence[str],
    pre_metabolites: Sequence[str],
    post_model: cobra.Model,
):
    """
    Run basic checks on new models (checking for objective flux).
    """

    # Check for objective flux
    new_genes = len(post_model.genes)
    new_rxn_ids = set([x.id for x in post_model.reactions]).difference(pre_reactions)
    new_cpd_ids = set([x.id for x in post_model.metabolites]).difference(
        pre_metabolites
    )
    test_flux = post_model.slim_optimize()

    # Report to user
    print(
        f"\tDraft reconstruction had {new_genes} genes, {len(pre_reactions)} reactions, and {len(pre_metabolites)} metabolites"
    )
    print(
        f"\tGapfilled {len(new_rxn_ids)} reactions and {len(new_cpd_ids)} metabolites\n"
    )
    print(
        f"\tFinal reconstruction has {len(post_model.reactions)} reactions and {len(post_model.metabolites)} metabolites"
    )
    print(f"\tFinal objective flux is {round(test_flux, 3)}")

    return test_flux
