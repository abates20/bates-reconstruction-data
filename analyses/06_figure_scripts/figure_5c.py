"""
figure_5c.py
------------

Code to generate the plot for Figure 5C.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "data",
    nargs="+",
    help="Paths to csvs with substrate utilization predictions"
)
parser.add_argument(
    "--exp-data",
    help="Path to csv with experimental substrate utilization data"
)
parser.add_argument(
    "--set-info",
    required=True,
    nargs="+",
    help="Path to csv with resource set info"
)
parser.add_argument(
    "--model-info",
    required=True,
    help="Path to model info csv"
)
parser.add_argument(
    "--out",
    default="figure_5c.svg",
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

from bates_reconstruction_data import palettes


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_5c_growth_metrics.csv"
    METRICS_DATA_EXISTS = METRICS_DATA.exists()
else:
    METRICS_DATA = None
    METRICS_DATA_EXISTS = False


def load_set_data(path: str) -> pl.DataFrame:
    """
    Load a csv with annotation set information.

    Assumes that the csv is located in a folder where the folder name is the
    organism_id.
    """
    org_id = Path(path).parent.name
    df = pl.read_csv(
        path,
        schema_overrides={"set_id": pl.String}
    )
    return df.with_columns(organism_id=pl.lit(org_id))


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
        schema_overrides={"set_id": pl.String}
    )

# Otherwise compute data for plotting
else:

    # Load experimental data
    experimental_data = (
        pl.scan_csv(args.exp_data)
        .select("organism_id", "compound", "growth")
    )

    # Calculate metrics for each input reaction data file
    metrics_pieces = []
    for path in args.data:

        # Load reaction data
        growth_data = pl.scan_csv(
            path,
            schema_overrides={"set_id": pl.String}
        )

        # Count number of substrates that support growth for each model
        growth_data = (
            growth_data
            .join(
                experimental_data,
                on=["organism_id", "compound"],
                suffix="_experimental"
            )
            .with_columns(
                label=(
                    pl.when(pl.col.growth & pl.col.growth_experimental)
                    .then(pl.lit("TP"))
                    .when(pl.col.growth)
                    .then(pl.lit("FP"))
                    .when(pl.col.growth_experimental)
                    .then(pl.lit("FN"))
                    .otherwise(pl.lit("TN"))
                )
            )
            .group_by("organism_id", "set_id", "label")
            .len()
            .pivot(
                "label",
                ["TP", "FP", "FN"],
                index=["organism_id", "set_id"],
                values="len"
            )
            .fill_null(0)
            .with_columns(
                precision = pl.col.TP / (pl.col.TP + pl.col.FP),
                recall = pl.col.TP / ( pl.col.TP + pl.col.FN)
            )
            .with_columns(
                f1_score = 2 * pl.col.precision * pl.col.recall / (pl.col.precision + pl.col.recall)
            )
            .drop("TP", "FP", "FN")
        )

        metrics_pieces.append(growth_data.collect())

    metrics = (
        pl.concat(metrics_pieces)
        .unpivot(
            on=["f1_score", "precision", "recall"],
            index=["organism_id", "set_id"],
            variable_name="metric",
            value_name="value"
        )
    )

    # Cache metrics if specified
    if args.cache is not None:
        metrics.write_csv(METRICS_DATA)


# Load info about resource sets and get sets with 0% random annotations
set_data = pl.concat([load_set_data(x) for x in args.set_info])
baseline_sets = set_data.filter(pl.col.gene.is_null())

# Prepare plot data by adding metadata
data = (
    metrics
    .join(set_data, on=["organism_id", "set_id"])
    .join(
        pl.read_csv(args.model_info),
        on="organism_id"
    )
    .join(
        metrics.join(baseline_sets, on=["organism_id", "set_id"], how="semi"),
        on=["organism_id", "metric"],
        suffix="_baseline"
    )
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}),
        pl.col.gene.str.split(".").list.last(),
        species_abbrev=species_abbreviation("species"),
        effect=pl.col.value - pl.col.value_baseline,
    )
)

# Grab the most impactful genes for each organism (biggest average effect on F1 score)
to_plot = (
    data
    .filter(pl.col.metric == "F1 Score")
    .group_by("organism_id", "gene")
    .agg(pl.col.effect.abs().mean())
    .group_by("organism_id")
    .agg(pl.col.gene, pl.col.effect, max_effect=pl.col.effect.max())
    .explode(["gene", "effect"], empty_as_null=False)
    .filter(pl.col.effect == pl.col.max_effect)
    .select("organism_id", "gene")
    .unique()
)

data = data.join(
    to_plot,
    on=["organism_id", "gene"],
    how="semi"
)

# Also set order of metrics to put F1 score on top
metric_order = ["Precision", "Recall", "F1 Score"]
data = (
    data.with_columns(
        pl.col.metric.cast(pl.Enum(metric_order))
    )
    .sort("metric")
)

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="gene",
            y="effect",
            color="metric"
        )
    )
    + pn.geom_hline(yintercept=0, linetype="dashed")
    + pn.geom_segment(
        pn.aes(xend="gene"),
        yend=0,
        position=pn.position_dodge(width=0.4),
        show_legend=False
    )
    + pn.geom_point(
        pn.aes(size="metric"),
        alpha=0.5,
        position=pn.position_dodge(width=0.4),
        show_legend=False
    )
    + pn.scale_color_manual(
        palettes.METRIC_COLOR,
        limits=["F1 Score", "Precision", "Recall"]
    )
    + pn.scale_size_manual(
        {"F1 Score": 2, "Precision": 1, "Recall": 1},
        limits=["F1 Score", "Precision", "Recall"]
    )
    + pn.labs(
        subtitle="Substrate utilization",
        x="",
        y="$\Delta$ metric score",
    )
    + pn.coord_flip()
    + pn.theme_classic()
    + pn.theme(
        figure_size=(2.9, 2.3),
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        axis_text_y=pn.element_text(size=8),
        plot_subtitle=pn.element_text(weight="bold", size=10),
        panel_spacing=0.005,
    )
)

plot.save(args.out)
