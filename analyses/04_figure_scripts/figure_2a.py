"""
figure_2a.py
------------

Code to generate the plot for Figure 2A.
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
    default="figure_2a.svg",
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
from sklearn.manifold import MDS

from bates_reconstruction_data import palettes


if args.cache is not None:
    MDS_DATA = Path(args.cache) / "figure_2a_MDS.csv"
    DATA_EXISTS = MDS_DATA.exists()
else:
    MDS_DATA = None
    DATA_EXISTS = None


if args.cache is not None and DATA_EXISTS:
    mds_data = pl.read_csv(MDS_DATA)        # ty:ignore[invalid-argument-type]

else:
    rxn_data = (
        pl.read_csv(args.data)
        .filter(pl.col.source_id != "BiGG")
    )
    X = (
        rxn_data
        .drop("organism_id", "source_id")
        .to_numpy()
    )

    # Run MDS
    n = 2
    method = MDS(
        n_components=n,
        metric_mds=False,       # non-metric for non-linear relationships
        metric="jaccard",       # Jaccard distance because all values are 0 or 1
        normalized_stress=True,
        init="classical_mds"
    )
    mds_data = (
        pl.DataFrame(
            method.fit_transform(X),
            schema={f"mds_dim{i+1}": pl.Float64 for i in range(n)}
        )
        .select(
            rxn_data.get_column("organism_id"),
            rxn_data.get_column("source_id"),
            pl.all()
        )
    )

    print(f"MDS normalized stress = {method.stress_}")

    # Save mds results
    if args.cache is not None:
        mds_data.write_csv(MDS_DATA)


# Combine MDS results with metadata
data = mds_data.join(
    pl.read_csv(args.set_info),
    on="source_id"
)

# Plot MDS reuslts
plot = (
    pn.ggplot(
        data,
        pn.aes(
            x="mds_dim1",
            y="mds_dim2",
            color="data_source",
            shape="reconstruction_tool"
        )
    )
    + pn.geom_point(alpha=0.7, size=4)
    + pn.scale_color_manual(palettes.DATA_SOURCE_COLOR)
    + pn.scale_shape_manual(palettes.RECONSTRUCTION_TOOL_SHAPE)
    + pn.labs(
        x="NMDS Dimension 1",
        y="NMDS Dimension 2",
        color="Data source",
        shape="Reconstruction tool"
    )
    + pn.coord_fixed(ratio=1)
    + pn.theme_classic(base_size=12)
    + pn.theme(
        figure_size=(3.5, 4.5),
        axis_ticks=pn.element_blank(),
        axis_text=pn.element_blank(),
        legend_position="bottom",
        legend_direction="vertical",
        legend_box_spacing=0.2,
        legend_spacing=25
    )
)

plot.save(args.out)
