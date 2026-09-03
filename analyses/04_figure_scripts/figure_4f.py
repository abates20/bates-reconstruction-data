"""
figure_4f.py
------------

Code to generate the plot for Figure 4F.
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
    default="figure_4f.svg",
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
from sklearn.metrics import precision_score, f1_score, recall_score
from sklearn.linear_model import LinearRegression


from bates_reconstruction_data import palettes, taxonomy


if args.cache is not None:
    REGRESSION_DATA = Path(args.cache) / "figure_4f_regression_results.csv"
    REGRESSION_DATA_EXISTS = REGRESSION_DATA.exists()

    DIST_DATA = Path(args.cache) / "figure_4f_taxonomic_distances.csv"
    DIST_DATA_EXISTS = DIST_DATA.exists()

    TAX_DMP = Path(args.cache) / "tax.dmp"
else:
    REGRESSION_DATA = None
    REGRESSION_DATA_EXISTS = False

    DIST_DATA = None
    DIST_DATA_EXISTS = False

    TAX_DMP = None


def linear_regression(data: dict[str, list[float]], x, y, scale=True) -> tuple[float]:
    """
    Run linear regression on a set of data and return the estimated
    coefficients.

    x is the name(s) of the independent variable(s). y is the name of the
    dependent variable.
    """
    def standard_scale(x):
        expr = pl.col(x)
        return (expr - expr.mean()) / expr.std()

    df = pl.DataFrame(data)

    # Scale all the data if scale is True
    if scale:
        df = df.select(
            standard_scale(x),
            standard_scale(y)
        )

    # Run linear regression
    X_ = df.select(x).to_numpy()
    y_ = df.get_column(y).to_numpy()
    model = LinearRegression(copy_X=False)
    model.fit(X_, y_)
    return tuple(model.coef_)


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


# Load computed regression results if cached
if args.cache is not None and REGRESSION_DATA_EXISTS:
    regression = pl.read_csv(REGRESSION_DATA)                              # ty:ignore[invalid-argument-type]

# Otherwise compute regression data for plotting
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

    x_vars = ["taxonomic_distance", "bigg_coverage"]
    y_var = "value"
    regression = (
        growth_metrics
        .join(tax_dists, on=["organism_id", "source_id"])
        .join(set_data, on="source_id")
        .group_by("organism_id", "metric")
        .agg(*x_vars, y_var)
        .with_columns(
            results=(
                pl.struct(*x_vars, y_var)
                .map_elements(
                    lambda d: linear_regression(d, x_vars, y_var),
                    return_dtype=pl.Struct({x: pl.Float64 for x in x_vars})
                )
            )
        )
        .drop(*x_vars, y_var)
        .unnest("results")
        .unpivot(
            x_vars,
            index=["organism_id", "metric"],
            variable_name="taxonomic_param",
            value_name="coeff"
        )
    )

    # Cache regression results if specified
    if args.cache is not None:
        regression.write_csv(REGRESSION_DATA)


# Prepare plot data by reformatting names
data = (
    regression
    .with_columns(
        pl.col.metric.replace(
            {"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"}
        ),
        pl.col.taxonomic_param.replace(
            {"taxonomic_distance": "Taxonomic\ndistance", "bigg_coverage": "Taxonomic\ncoverage"}
        )
    )
)

DODGE = 0.7
plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="taxonomic_param",
            y="coeff",
            color="metric"
        )
    )
    + pn.geom_hline(
        yintercept=0,
        linetype="dashed"
    )
    + pn.geom_point(
        alpha=0.2,
        size=1,
        position=pn.position_jitterdodge(
            jitter_width=0.1,
            dodge_width=DODGE
        ),
        show_legend=False
    )
    + pn.stat_summary(
        geom=pn.geom_errorbar,
        width=0.2,
        size=0.7,
        fun_data="mean_cl_boot",
        fun_args=dict(confidence_interval=0.95),
        position=pn.position_dodge(width=DODGE)
    )
    + pn.stat_summary(
        geom=pn.geom_errorbar,
        width=0.5,
        size=1.2,
        fun_ymin=np.mean,
        fun_ymax=np.mean,
        position=pn.position_dodge(width=DODGE)
    )
    + pn.scale_color_manual(palettes.METRIC_COLOR)
    + pn.labs(
        subtitle="Substrate utilization",
        x="",
        y="Effect on\nprecision/recall",
        color=""
    )
    + pn.theme_classic()
    + pn.theme(
        figure_size=(3,2.2),
        plot_subtitle=pn.element_text(weight="bold", size=10),
    )
)

plot.save(args.out)
