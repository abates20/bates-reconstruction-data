"""
figure_3a.py
------------

Code to generate the plot for Figure 3A.
"""


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
    default="figure_3a.svg",
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
    .select("database", "species")
)
bigg_data = (
    pl.read_csv(args.bigg_data)
    .with_columns(database=pl.lit("BiGG"))
    .select("database", "species")
)

# Combine and get number of species in each db
data = (
    pl.concat([kegg_data, bigg_data])
    .group_by("database")
    .agg(pl.col.species.unique().len())
)

# Plot species counts
plot = (
    pn.ggplot(
        data,
        pn.aes(x="database", y="species", fill="database")
    )
    + pn.geom_col(width=0.8, show_legend=False)
    + pn.geom_text(
        pn.aes(label="species"),
        nudge_y=200,
        ha="left",
        size=8
    )
    + pn.scale_fill_manual(palettes.DATABASE_COLOR)
    + pn.labs(
        x="",
        y="Number of bacteria\nspecies"
    )
    + pn.coord_flip(ylim=(0, 11000))
    + pn.theme_classic()
    + pn.theme(figure_size=(2.1,1.4))
)

plot.save(args.out)
