"""
figure_5a.py
------------

Code to generate the plot for Figure 5A.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "data",
    nargs="+",
    help="Paths to csvs with reaction data"
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
    default="figure_5a.svg",
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
    METRICS_DATA = Path(args.cache) / "figure_5a_reaction_metrics.csv"
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
    # Calculate metrics for each input reaction data file
    metrics_pieces = []
    for path in args.data:

        # Load reaction data
        rxn_data = pl.scan_csv(
            path,
            schema_overrides={"set_id": pl.String}
        )

        # Pull out the data for the BiGG model
        bigg_data = (
            rxn_data
            .filter(pl.col.set_id == "BiGG")
            .unpivot(
                index=["organism_id", "set_id"],
                variable_name="reaction",
                value_name="is_present"
            )
            .drop("organism_id", "set_id")
        )

        # Calculate metrics using BiGG model as ground truth
        rxn_data = (
            rxn_data
            .filter(pl.col.set_id != "BiGG")
            .unpivot(
                index=["organism_id", "set_id"],
                variable_name="reaction",
                value_name="is_present"
            )
            .join(bigg_data, on="reaction", suffix="_BiGG")
            .with_columns(
                label=(
                    pl.when(pl.col.is_present.eq(1) & pl.col.is_present_BiGG.eq(1))
                    .then(pl.lit("TP"))
                    .when(pl.col.is_present.eq(1))
                    .then(pl.lit("FP"))
                    .when(pl.col.is_present_BiGG.eq(1))
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

        metrics_pieces.append(rxn_data.collect())

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
        species_abbrev=species_abbreviation("species")
    )
    .with_columns(
        close_to_baseline=(pl.col.value - pl.col.value_baseline).abs().lt(0.01)     # 0.01 as cutoff to determine what points are "close to the baseline"
    )
)

# Downsample the points that are close to the baseline value (to reduce over plotting)
to_plot = (
    data
    .filter(pl.col.close_to_baseline)
    .group_by("organism_id", "metric")
    .agg(pl.col.set_id.sort().gather_every(20))
    .explode("set_id", empty_as_null=False)
    .vstack(
        data.filter(~pl.col.close_to_baseline).select("organism_id", "metric", "set_id")
    )
)

data = data.join(
    to_plot,
    on=["organism_id", "metric", "set_id"],
    how="semi"
)

# Set order of species so that they line up with figure 5b
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
data = (
    data.with_columns(
        pl.col.species_abbrev.cast(pl.Enum(species_order)),
        pl.col.metric.cast(pl.Enum(metric_order)),
        pl.col.set_id.cast(pl.Int64)    # Puts set_ids in order when plotting
    )
    .sort("metric")
)

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="set_id",
            y="value",
            color="metric"
        )
    )
    + pn.facet_wrap("species_abbrev", ncol=5, scales="free_x")
    + pn.geom_point(
        alpha=0.5,
        size=0.5,
    )
    + pn.scale_color_manual(
        palettes.METRIC_COLOR,
        limits=["F1 Score", "Precision", "Recall"],
        guide=None
    )
    + pn.labs(
        subtitle="Reaction content",
        x="",
        y="Metric score",
        color="",
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic()
    + pn.theme(
        figure_size=(4, 2.5),
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        axis_text_x=pn.element_blank(),
        axis_ticks_x=pn.element_blank(),
        legend_title=pn.element_blank(),
        legend_text=pn.element_text(size=11),
        legend_key_spacing_y=10,
        plot_subtitle=pn.element_text(weight="bold", size=10),
        panel_spacing=0.005
    )
)

plot.save(args.out)
