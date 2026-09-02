"""
figure_4c.py
------------

Code to generate the plot for Figure 4C.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "data",
    help="Path to csv with reaction data"
)
parser.add_argument(
    "--set-info",
    required=True,
    help="Path to csv with resource set info"
)
parser.add_argument(
    "--model-info",
    required=True,
    help="Path to model info csv"
)
parser.add_argument(
    "--out",
    default="figure_4c.svg",
    help="Path to save figure"
)
parser.add_argument(
    "--cache",
    help="Path to directory for caching computed files"
)
args = parser.parse_args()


# Additional imports
from pathlib import Path
from collections import defaultdict

import polars as pl
import plotnine as pn
from sklearn.metrics import precision_score, f1_score, recall_score

from bates_reconstruction_data import palettes, taxonomy


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_4c_reaction_metrics.csv"
    METRICS_DATA_EXISTS = METRICS_DATA.exists()

    DIST_DATA = Path(args.cache) / "figure_4c_taxonomic_distances.csv"
    DIST_DATA_EXISTS = DIST_DATA.exists()

    TAX_DMP = Path(args.cache) / "tax.dmp"
else:
    METRICS_DATA = None
    METRICS_DATA_EXISTS = False

    DIST_DATA = None
    DIST_DATA_EXISTS = False

    TAX_DMP = None


def species_abbreviation(name: str) -> pl.Expr:
    """
    Polars expression to abbreviate the name of a species (i.e., first initial
    of the genus followed by species name).
    """
    name_pieces = (
        pl.col(name)
        .str.splitn(" ", 2)
        .struct.rename_fields(["genus", "species"])
    )
    genus = name_pieces.struct.field("genus")
    species = name_pieces.struct.field("species")
    return genus.str.head(1) + ". " + species


# Load computed metrics if cached
if args.cache is not None and METRICS_DATA_EXISTS:
    metrics = pl.read_csv(
        METRICS_DATA,                                       # ty:ignore[invalid-argument-type]
        schema_overrides={"source_id": pl.String}
    )

# Otherwise compute data for plotting
else:
    # Load reaction data
    rxn_data = pl.read_csv(
        args.data,
        schema_overrides={"source_id": pl.String}
    )

    # Index the reaction data
    data_index = defaultdict(dict)
    row_iter = rxn_data.select("organism_id", "source_id").iter_rows()
    for i, (org_id, source_id) in enumerate(row_iter):
        data_index[org_id][source_id] = i

    # Get source/organism combinations
    organisms_and_sources = (
        rxn_data
        .filter(pl.col.source_id != "BiGG")
        .select("organism_id", "source_id")
        .unique()
        .iter_rows()
    )

    # Dictionary to compile metrics
    metric_dict = {
        "organism_id": [],
        "source_id": [],
        "precision": [],
        "recall": [],
        "f1_score": []
    }

    # Calculate metrics
    matrix = rxn_data.drop("organism_id", "source_id")
    for org_id, source_id in organisms_and_sources:
        pred = matrix.row(data_index[org_id][source_id])
        true = matrix.row(data_index[org_id]["BiGG"])

        metric_dict["organism_id"].append(org_id)
        metric_dict["source_id"].append(source_id)
        metric_dict["precision"].append(precision_score(true, pred))
        metric_dict["recall"].append(recall_score(true, pred))
        metric_dict["f1_score"].append(f1_score(true, pred))

    metrics = (
        pl.DataFrame(metric_dict)
        .unpivot(
            on=["f1_score", "precision", "recall"],
            index=["organism_id", "source_id"],
            variable_name="metric",
            value_name="value"
        )
    )

    # Cache metrics if specified
    if args.cache is not None:
        metrics.write_csv(METRICS_DATA)


# Load info about resource sets
set_data = pl.read_csv(
    args.set_info,
    schema_overrides={"source_id": pl.String}
)

# Load model info
model_info = pl.read_csv(args.model_info)


# Load computed taxonomic distances if cached
if args.cache is not None and DIST_DATA_EXISTS:
    tax_dists = pl.read_csv(
        DIST_DATA,                                       # ty:ignore[invalid-argument-type]
        schema_overrides={"source_id": pl.String}
    )

# Otherwise compute taxonomic distances
else:
    model_taxa = taxonomy.get_taxa(model_info.get_column("species"), TAX_DMP)

    dist_dict = {
        "source_id": [],
        "organism_id": [],
        "taxonomic_distance": []
    }

    set_iter = (
        set_data
        .select(
            "source_id",
            pl.col.selected_species.str.split(",")
        )
        .iter_rows()
    )

    for source_id, species_list in set_iter:
        set_taxa = taxonomy.get_taxa(species_list, TAX_DMP)

        for org_id, model_tx in zip(model_info.get_column("organism_id"), model_taxa):
            dist = min(taxonomy.taxonomic_distance(model_tx, tx) for tx in set_taxa)
            dist_dict["source_id"].append(source_id)
            dist_dict["organism_id"].append(org_id)
            dist_dict["taxonomic_distance"].append(dist)

    tax_dists = pl.DataFrame(dist_dict)

    # Cache distances if specified
    if args.cache is not None:
        tax_dists.write_csv(DIST_DATA)


# Prepare plot data by adding metadata
data = (
    metrics
    .join(set_data, on="source_id")
    .join(model_info, on="organism_id")
    .join(tax_dists, on=["source_id", "organism_id"])
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}),
        species_abbrev=species_abbreviation("species")
    )
)

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="taxonomic_distance",
            y="value",
            color="metric"
        )
    )
    + pn.facet_wrap("species_abbrev", nrow=2)
    + pn.geom_point(
        alpha=0.3,
        size=0.7
    )
    + pn.geom_path(
        pn.aes(group="metric", linetype="metric", size="metric"),
        stat="summary"
    )
    + pn.scale_color_manual(palettes.METRIC_COLOR, guide=None)
    + pn.scale_linetype_manual(
        {"F1 Score": "solid", "Precision": "dashed", "Recall": "dashed"},
        guide=None
    )
    + pn.scale_size_manual(
        {"F1 Score": 1.25, "Precision": 0.6, "Recall": 0.6},
        guide=None
    )
    + pn.labs(
        subtitle="Reactions",
        x="Taxonomic distance",
        y="",
        linetype="",
        size=""
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic()
    + pn.theme(
        figure_size=(3.7,2.4),
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        plot_subtitle=pn.element_text(weight="bold", size=10),
        panel_spacing=0.01
    )
)

plot.save(args.out)
