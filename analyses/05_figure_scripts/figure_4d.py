"""
figure_4d.py
------------

Code to generate the plot for Figure 4D.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "data",
    help="Path to csv with substrate utilization data/predictions"
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
    default="figure_4d.svg",
    help="Path to save figure"
)
parser.add_argument(
    "--cache",
    help="Path to directory for caching computed files"
)
args = parser.parse_args()


# Additional imports
from pathlib import Path

import polars as pl
import plotnine as pn
import numpy as np
from scipy import stats
from sklearn.metrics import precision_score, f1_score, recall_score
from mizani.labels import label_percent

from bates_reconstruction_data import palettes


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_4d_growth_metrics.csv"
    METRICS_DATA_EXISTS = METRICS_DATA.exists()

    STATS_DATA = Path(args.cache) / "figure_4d_metric_stats.csv"
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
    growth_metrics = pl.read_csv(
        METRICS_DATA,                                       # ty:ignore[invalid-argument-type]
        schema_overrides={"source_id": pl.String}
    )

# Otherwise compute data for plotting
else:
    # Load growth data
    growth_data = (
        pl.read_csv(
            args.data,
            schema_overrides={"source_id": pl.String}
        )
        .select("organism_id", "source_id", "compound", "growth")
    )

    # Create separate df per organism (to remove compounds missing experimental data)
    growth_dfs: dict[str, pl.DataFrame] = {}
    for org_id in growth_data.get_column("organism_id").unique():
        df = (
            growth_data
            .filter(pl.col.organism_id == org_id)
            .pivot("source_id", values="growth")
            .drop("organism_id", "compound")
            .drop_nulls()
        )
        growth_dfs[org_id] = df

    # Get source/organism combinations
    organisms_and_sources = (
        growth_data
        .filter(pl.col.source_id != "experimental")
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
    for org_id, source_id in organisms_and_sources:
        df = growth_dfs[org_id]
        pred = df.get_column(source_id)
        true = df.get_column("experimental")

        metric_dict["organism_id"].append(org_id)
        metric_dict["source_id"].append(source_id)
        metric_dict["precision"].append(precision_score(true, pred))
        metric_dict["recall"].append(recall_score(true, pred))
        metric_dict["f1_score"].append(f1_score(true, pred))

    growth_metrics = (
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
        growth_metrics.write_csv(METRICS_DATA)


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
        growth_metrics
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
    metric_stats
    .join(
        pl.read_csv(args.model_info),
        on="organism_id"
    )
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}),
        species_abbrev=species_abbreviation("species")
    )
)

# Set order of metrics to put F1 score first
metric_order = ["F1 Score", "Precision", "Recall"]
data = data.with_columns(pl.col.metric.cast(pl.Enum(metric_order)))

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="metric",
            y="slope",
            color="metric",
        )
    )
    + pn.geom_hline(
        yintercept=0,
        linetype="dashed"
    )
    + pn.geom_point(
        alpha=0.2,
        position=pn.position_jitter(width=0.1, height=0)
    )
    + pn.stat_summary(
        geom=pn.geom_errorbar,
        width=0.15,
        size=1,
        fun_data="mean_cl_boot",
        fun_args=dict(confidence_interval=0.95)
    )
    + pn.stat_summary(
        geom=pn.geom_errorbar,
        width=0.4,
        size=1.5,
        fun_ymin=np.mean,
        fun_ymax=np.mean
    )
    + pn.scale_color_manual(palettes.METRIC_COLOR, guide=None)
    + pn.labs(
        subtitle="Substrate utilization",
        x="",
        y="Regression coefficient",
        color="",
        linetype="",
        size=""
    )
    + pn.theme_classic()
    + pn.theme(
        figure_size=(1.8, 1.7),
        axis_title_y=pn.element_text(size=10),
        axis_text_x=pn.element_text(rotation=45, size=8, margin={"t": 0, "unit": "pt"}),
        plot_subtitle=pn.element_text(weight="bold", size=10),
    )
)

plot.save(args.out)
