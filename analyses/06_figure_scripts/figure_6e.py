"""
figure_6e.py
------------

Code to generate the plot for Figure 6E.
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
    default="figure_6e.svg",
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
from sklearn.metrics import precision_score, f1_score, recall_score

from bates_reconstruction_data import palettes, taxonomy


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_6e_growth_metrics.csv"
    METRICS_DATA_EXISTS = METRICS_DATA.exists()

    DIST_DATA = Path(args.cache) / "figure_6e_taxonomic_distances.csv"
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
    growth_metrics
    .join(set_data, on="source_id")
    .join(model_info, on="organism_id")
    .join(tax_dists, on=["source_id", "organism_id"])
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}),
        species_abbrev=species_abbreviation("species")
    )
)

# Calculate spacing values as fraction of total width (to make spacing same as 6c)
width = 2.5
height = 2.3
legend_box_spacing = 0.04 / width
panel_spacing = 0.02 / width
plot_margin = 0.025 / width

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="taxonomic_distance",
            y="value",
            color="metric",
        )
    )
    + pn.facet_wrap("species_abbrev", nrow=2)
    + pn.geom_point(
        alpha=0.3,
        size=0.7,
    )
    + pn.geom_path(
        pn.aes(group="metric", linetype="metric", size="metric"),
        stat="summary"
    )
    + pn.scale_color_manual(palettes.METRIC_COLOR)
    + pn.scale_linetype_manual({"F1 Score": "solid", "Precision": "dashed", "Recall": "dashed"})
    + pn.scale_size_manual({"F1 Score": 1.25, "Precision": 0.6, "Recall": 0.6})
    + pn.labs(
        subtitle="Substrate utilization",
        x="Taxonomic distance",
        y="Metric score",
        color="",
        linetype="",
        size=""
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic()
    + pn.theme(
        figure_size=(width, height),
        plot_margin=plot_margin,
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        axis_title_x=pn.element_text(margin={"t": 5 / width}),
        legend_title=pn.element_blank(),
        legend_position="bottom",
        legend_justification="right",
        legend_box_spacing=legend_box_spacing,
        plot_subtitle=pn.element_text(weight="bold", size=10, margin={"t": 0, "b": 5 / width}),
        panel_spacing=panel_spacing,
    )
)

plot.save(args.out)
