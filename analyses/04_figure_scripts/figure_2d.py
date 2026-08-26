"""
figure_2d.py
------------

Code to generate the plot for Figure 2D.
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
    "--out",
    default="figure_2d.svg",
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
from sklearn.metrics import precision_score, f1_score, recall_score

from bates_reconstruction_data import palettes


if args.cache is not None:
    METRICS_DATA = Path(args.cache) / "figure_2d_growth_metrics.csv"
    DATA_EXISTS = METRICS_DATA.exists()
else:
    METRICS_DATA = None
    DATA_EXISTS = None


if args.cache is not None and DATA_EXISTS:
    growth_metrics = pl.read_csv(METRICS_DATA)

else:
    # Load growth data
    growth_data = (
        pl.read_csv(args.data)
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

    growth_metrics = pl.DataFrame(metric_dict)

    # Cache metrics if specified
    if args.cache is not None:
        growth_metrics.write_csv(METRICS_DATA)


# Prepare plot data by pivoting to long format and adding metadata
data = (
    growth_metrics
    .unpivot(
        on=["f1_score", "precision", "recall"],
        index=["organism_id", "source_id"],
        variable_name="metric",
        value_name="value"
    )
    .with_columns(
        pl.col.metric.replace({"precision": "Precision", "recall": "Recall", "f1_score": "F1 Score"})
    )
    .join(
        pl.read_csv(args.set_info),
        on="source_id"
    )
)

plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="reconstruction_tool",
            y="value",
            color="data_source",
            shape="reconstruction_tool"
        )
    )
    + pn.facet_wrap("metric")
    + pn.geom_point(
        alpha=0.15,
        size=0.8,
        position=pn.position_jitterdodge(dodge_width=0.7, jitter_width=0.1),
        show_legend=False
    )
    + pn.geom_errorbar(
        width=0.2,
        size=0.8,
        stat="summary",
        position=pn.position_dodge(width=0.7),
        show_legend=False
    )
    + pn.geom_errorbar(
        width=0.5,
        size=1.2,
        stat="summary",
        fun_ymin=np.mean,
        fun_ymax=np.mean,
        position=pn.position_dodge(width=0.7),
        show_legend=False
    )
    + pn.scale_color_manual(palettes.DATA_SOURCE_COLOR)
    + pn.scale_shape_manual(palettes.RECONSTRUCTION_TOOL_SHAPE)
    + pn.labs(
        x="",
        y="",
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic()
    + pn.theme(
        figure_size=(3,2.3),
        axis_text_x=pn.element_text(rotation="vertical")
    )
)

plot.save(args.out)
