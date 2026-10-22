"""
figure_4b.py
------------

Code to generate the plot for Figure 4B.
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
    default="figure_4b.svg",
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
import numpy as np
from scipy import stats
from sklearn.metrics import precision_score, f1_score, recall_score

from bates_reconstruction_data import palettes


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_4b_reaction_metrics.csv"
    METRICS_DATA_EXISTS = METRICS_DATA.exists()

    STATS_DATA = Path(args.cache) / "figure_4b_metric_stats.csv"
    STATS_DATA_EXISTS = STATS_DATA.exists()
else:
    METRICS_DATA = None
    METRICS_DATA_EXISTS = False

    STATS_DATA = None
    STATS_DATA_EXISTS = False


def linear_regression(data: dict):
    """
    Run linear regression.

    `data` is a dictionary with keys 'x' (the independent variable) and 'y' (the
    dependent variable).

    Returns the slope and pvalue.
    """
    res = stats.linregress(data["x"], data["y"])
    slope = None if np.isnan(res.slope) else res.slope
    pval = None if np.isnan(res.pvalue) else res.pvalue
    return slope, pval


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


# Load info about resource sets and get sets with 0% random annotations
set_data = pl.read_csv(
    args.set_info,
    schema_overrides={"source_id": pl.String}
)
baseline_sets = set_data.filter(pl.col.fraction == 0)


# Load computed stats if cached
if args.cache is not None and STATS_DATA_EXISTS:
    metric_stats = pl.read_csv(STATS_DATA)      # ty:ignore[invalid-argument-type]

# Otherwise run stats
else:

    # Run linear regression for each combo of species and metric
    metric_stats = (
        metrics
        .join(set_data, on="source_id")
        .group_by("organism_id", "metric")
        .agg("fraction", "value")
        .with_columns(
            result=(
                pl.struct(x="fraction", y="value")
                .map_elements(
                    linear_regression,
                    return_dtype=pl.Struct({"slope": pl.Float64, "pval": pl.Float64})
                )
            )
        )
        .unnest("result")
        .drop("fraction", "value")
    )

    # False discovery control to get adjusted p-values
    adj_pvals = (
        metric_stats
        .drop_nulls()
        .select(
            "organism_id",
            "metric",
            adj_pval=pl.col.pval.map_batches(
                stats.false_discovery_control,
                return_dtype=pl.Float64
            )
        )
    )
    metric_stats = metric_stats.join(
        adj_pvals,
        on=["organism_id", "metric"],
        how="left"
    )

    # Cache stats if specified
    if args.cache is not None:
        metric_stats.write_csv(STATS_DATA)


# Prepare plot data by adding metadata
data = (
    metrics
    .join(set_data, on="source_id")
    .join(
        pl.read_csv(args.model_info),
        on="organism_id"
    )
    .join(
        metric_stats,
        on=["organism_id", "metric"],
        how="left"
    )
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}),
        species_abbrev=species_abbreviation("species")
    )
)

# Set order of species so that they line up with figure 4c
# Also set order of metrics to put F1 score on top
species_order = [
    "A. baumannii",
    "B. subtilis",
    "C. difficile",
    "C. ljungdahlii",
    "L. cremoris",
    "E. coli",
    "P. putida",
    "S. aureus",
    "S. sonnei",
    "Y. pestis"
]
metric_order = ["Precision", "Recall", "F1 Score"]
data = data.with_columns(
    pl.col.species_abbrev.cast(pl.Enum(species_order)),
    pl.col.metric.cast(pl.Enum(metric_order))
)

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="fraction",
            y="value",
            color="metric"
        )
    )
    + pn.facet_wrap("species_abbrev", ncol=5)
    + pn.geom_point(
        alpha=0.1,
        size=0.7
    )
    + pn.geom_smooth(
        pn.aes(group="metric", linetype="metric", size="metric"),
        method="lowess",
        se=False
    )
    + pn.scale_color_manual(
        palettes.METRIC_COLOR,
        limits=["F1 Score", "Precision", "Recall"]
    )
    + pn.scale_linetype_manual(
        {"F1 Score": "solid", "Precision": "dashed", "Recall": "dashed"},
        limits=["F1 Score", "Precision", "Recall"]
    )
    + pn.scale_size_manual(
        {"F1 Score": 1.25, "Precision": 0.6, "Recall": 0.6},
        limits=["F1 Score", "Precision", "Recall"]
    )
    + pn.scale_x_continuous(breaks=[0.0, 0.3, 0.6, 0.9])
    + pn.labs(
        subtitle="Reaction content",
        x="Fraction of genes with random annotation",
        y="Metric score",
        color="",
        linetype="",
        size=""
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic()
    + pn.theme(
        figure_size=(4.1,3.2),
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        axis_text_x=pn.element_text(rotation=45, size=7, margin={"t": 0, "unit": "pt"}),
        legend_position="bottom",
        legend_title=pn.element_blank(),
        legend_text=pn.element_text(size=11),
        legend_key_spacing_x=20,
        plot_subtitle=pn.element_text(weight="bold", size=10),
        panel_spacing=0.005
    )
)

plot.save(args.out)
