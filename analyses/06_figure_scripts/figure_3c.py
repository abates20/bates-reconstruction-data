"""
figure_3c.py
------------

Code to generate the plot for Figure 3C.
"""
from plotnine.themes.themeable import legend_title_align


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "--kegg-data",
    help="Path to csv with KEGG data"
)
parser.add_argument(
    "--bigg-data",
    help="Path to csv with BiGG data"
)
parser.add_argument(
    "--out",
    default="figure_3c.svg",
    help="Path to save figure"
)
args = parser.parse_args()


# Additional imports
import polars as pl
import plotnine as pn

from bates_reconstruction_data import palettes


# Load data
kegg_data = (
    pl.read_csv(args.kegg_data)
    .with_columns(database=pl.lit("KEGG"))
    .select("database", "n_genes")
)
bigg_data = (
    pl.read_csv(args.bigg_data)
    .with_columns(database=pl.lit("BiGG"))
    .select("database", "n_genes")
)

# Combine data for each db
data = pl.concat([kegg_data, bigg_data])

# Plot species counts
plot = (
    pn.ggplot(
        data,
        pn.aes(x="n_genes", fill="database")
    )
    + pn.geom_density(alpha=0.8)
    + pn.scale_fill_manual(palettes.DATABASE_COLOR, limits=["KEGG", "BiGG"])
    + pn.labs(
        x="Annotated genes\nper organism",
        y="Density"
    )
    + pn.theme_classic()
    + pn.theme(
        figure_size=(2.2, 1.4),
        axis_text_y=pn.element_blank(),
        axis_ticks_y=pn.element_blank(),
        legend_title=pn.element_blank()
    )
)

plot.save(args.out)
