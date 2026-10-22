"""
figure_3b.py
------------

Code to generate the plot for Figure 3B.
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
    default="figure_3b.svg",
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
    .select("database", "reaction")
)
bigg_data = (
    pl.read_csv(args.bigg_data)
    .with_columns(database=pl.lit("BiGG"))
    .select("database", "reaction")
)

# Combine and get number of reactions in each db
data = (
    pl.concat([kegg_data, bigg_data])
    .group_by("database")
    .agg(pl.col.reaction.unique().len())
)

# Plot species counts
plot = (
    pn.ggplot(
        data,
        pn.aes(x="database", y="reaction", fill="database")
    )
    + pn.geom_col(width=0.8, show_legend=False)
    + pn.geom_text(
        pn.aes(label="reaction"),
        nudge_y=200,
        ha="left",
        size=8
    )
    + pn.scale_fill_manual(palettes.DATABASE_COLOR)
    + pn.labs(
        x="",
        y="Number of reactions\nlinked to genes"
    )
    + pn.coord_flip(ylim=(0, 9000))
    + pn.theme_classic()
    + pn.theme(figure_size=(2, 1.4))
)

plot.save(args.out)
