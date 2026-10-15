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

# Set width of figure elements to ensure panels align with figure 3b
total_width = 3.7 # in
buffer = 0.5 / 72 # in (when svg saved, there is a 0.5pt buffer added on each edge?)
panel_width = 0.6 # in
legend_width = 0.958 # in
panel_spacing = 0.01 * total_width # in
axis_text_y_width = 34.5 / 72 # in
plot_margin = 0.01 * total_width # in
legend_margin = 0.1 # in
nrow = 2 # number of rows of panels
n_panels = data.get_column("species_abbrev").unique().len() / nrow
legend_box_spacing = (
    total_width
    + (2 * buffer)
    - (2 * plot_margin)
    - (2 * legend_margin)
    - legend_width
    - (panel_width * n_panels)
    - (panel_spacing * (n_panels - 1))
    - axis_text_y_width
)

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="taxonomic_distance",
            y="value",
            color="metric",
        )
    )
    + pn.facet_wrap("species_abbrev", nrow=nrow)
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
        figure_size=(total_width,2.4),
        strip_text_x=pn.element_text(style="italic", family="arial", size=7),
        legend_title=pn.element_blank(),
        plot_subtitle=pn.element_text(weight="bold", size=10),
        panel_spacing=panel_spacing / total_width,                                  # convert inches to fraction
        legend_box_spacing=legend_box_spacing / total_width,                        # convert inches to fraction (plotnine docs say pts but seems to actually be fraction)
        legend_text=pn.element_text(size=12),
        legend_key_spacing_y=10,
        legend_box_margin=legend_margin * 72,                                       # convert inches to pts
        plot_margin_left=plot_margin / total_width,                                 # convert inches to fraction
        plot_margin_right=plot_margin / total_width,                                # convert inches to fraction
    )
)

plot.save(args.out)
