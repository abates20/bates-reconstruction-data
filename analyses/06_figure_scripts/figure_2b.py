"""
figure_2b.py
------------

Code to generate the plot for Figure 2B.
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
    "--out",
    default="figure_2b.svg",
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
    METRICS_DATA = Path(args.cache) / "figure_2b_reaction_metrics.csv"
    DATA_EXISTS = METRICS_DATA.exists()
else:
    METRICS_DATA = None
    DATA_EXISTS = None


def get_rxn_array(df: pl.DataFrame, org_id: str, source_id: str) -> tuple:
    row = (
        df.filter(
            pl.col.organism_id == org_id,
            pl.col.source_id == source_id
        )
        .drop("organism_id", "source_id")
        .row(0)
    )
    return row


if args.cache is not None and DATA_EXISTS:
    rxn_metrics = pl.read_csv(METRICS_DATA)     # ty:ignore[invalid-argument-type]

else:
    # Load reaction data
    rxn_data = pl.read_csv(args.data)

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
    for org_id, source_id in organisms_and_sources:
        pred = get_rxn_array(rxn_data, org_id, source_id)
        true = get_rxn_array(rxn_data, org_id, "BiGG")

        metric_dict["organism_id"].append(org_id)
        metric_dict["source_id"].append(source_id)
        metric_dict["precision"].append(precision_score(true, pred))
        metric_dict["recall"].append(recall_score(true, pred))
        metric_dict["f1_score"].append(f1_score(true, pred))

    rxn_metrics = pl.DataFrame(metric_dict)

    # Cache metrics if specified
    if args.cache is not None:
        rxn_metrics.write_csv(METRICS_DATA)


# Prepare plot data by pivoting to long format and adding metadata
data = (
    rxn_metrics
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
        subtitle="Reaction content",
        x="",
        y="Metric score",
    )
    + pn.coord_cartesian(ylim=(0, 1))
    + pn.theme_classic(base_family="Helvetica")
    + pn.theme(
        figure_size=(3,2),
        axis_text_x=pn.element_text(rotation=40),
        plot_subtitle=pn.element_text(weight="bold", size=10)
    )
)

# Hack to put the x-axis labels in the right place
# See: https://github.com/has2k1/plotnine/issues/1141#issuecomment-6347184892
fig = plot.draw()
fig.draw_without_rendering()
fig.set_layout_engine("none")
for ax in fig.axes:
    for label in ax.get_xticklabels():
        label.set(ha="right", va="top", rotation_mode="anchor", y=0)
fig.savefig(args.out)
