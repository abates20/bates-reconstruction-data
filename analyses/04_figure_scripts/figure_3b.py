"""
figure_3b.py
------------

Code to generate the plot for Figure 3B.
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
    default="figure_3b.svg",
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
from mizani.labels import label_percent

from bates_reconstruction_data import palettes


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_3b_reaction_metrics.csv"
    METRICS_DATA_EXISTS = METRICS_DATA.exists()

    STATS_DATA = Path(args.cache) / "figure_3b_metric_stats.csv"
    STATS_DATA_EXISTS = STATS_DATA.exists()
else:
    METRICS_DATA = None
    METRICS_DATA_EXISTS = False

    STATS_DATA = None
    STATS_DATA_EXISTS = False


def t_test(data: dict):
    """
    Run a 1 sample t-test.

    `data` is a dictionary with keys 'x' (the sample values) and 'mu' (the
    population mean that is being tested against).

    Returns the test statistic and pvalue.
    """
    res = stats.ttest_1samp(data["x"], data["mu"])
    stat = None if np.isnan(res.statistic) else res.statistic
    pval = None if np.isnan(res.pvalue) else res.pvalue
    return stat, pval


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

    # Get baseline metrics (sets with 0% random annotations)
    baseline_metrics = (
        metrics
        .join(baseline_sets, on="source_id", how="semi")
        .drop("source_id")
    )

    # Test sets with random annotations to see if metrics are different from baseline
    metric_stats = (
        metrics
        .join(baseline_sets, on="source_id", how="anti")
        .join(set_data, on="source_id")
        .group_by("organism_id", "metric", "fraction")
        .agg("value")
        .join(baseline_metrics, on=["organism_id", "metric"], suffix="_baseline")
        .with_columns(
            ttest_result=(
                pl.struct(x="value", mu="value_baseline")
                .map_elements(
                    t_test,
                    return_dtype=pl.Struct({"t_stat": pl.Float64, "pval": pl.Float64})
                )
            )
        )
        .unnest("ttest_result")
        .drop("value", "value_baseline")
    )

    # False discovery control to get adjusted p-values
    adj_pvals = (
        metric_stats
        .drop_nulls()
        .select(
            "organism_id",
            "metric",
            "fraction",
            adj_pval=pl.col.pval.map_batches(
                stats.false_discovery_control,
                return_dtype=pl.Float64
            )
        )
    )
    metric_stats = metric_stats.join(
        adj_pvals,
        on=["organism_id", "metric", "fraction"],
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
        on=["organism_id", "metric", "fraction"],
        how="left"
    )
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}),
        species_abbrev=species_abbreviation("species")
    )
)

# Get sorted fraction order
fraction_levels = (
    data
    .get_column("fraction")
    .unique()
    .sort()
)

# Split up data for 0% random annotations and >0% random annotations
baseline_data = data.join(baseline_sets, on="source_id", how="semi")
data = data.join(baseline_sets, on="source_id", how="anti")

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="factor(fraction, categories=fraction_levels)",
            y="value",
            color="metric"
        )
    )
    + pn.facet_grid(cols="species_abbrev")
    + pn.geom_point(data=baseline_data, shape="*", size=2)
    + pn.geom_point(
        alpha=0.3,
        size=0.7,
        position=pn.position_jitter(width=0.1)
    )
    + pn.geom_path(
        pn.aes(group="metric", linetype="metric", size="metric"),
        pl.concat([data, baseline_data]),
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
    + pn.scale_x_discrete(labels=label_percent())
    + pn.labs(
        subtitle="Reactions",
        x="Fraction of genes with random annotation",
        y="",
        linetype="",
        size=""
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic()
    + pn.theme(
        figure_size=(7,2.2),
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        axis_text_x=pn.element_text(rotation="vertical"),
        plot_subtitle=pn.element_text(weight="bold", size=10),
        panel_spacing=0.005
    )
)

plot.save(args.out)
