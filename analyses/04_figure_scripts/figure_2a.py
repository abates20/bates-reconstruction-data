"""
figure_2a.py
------------

Code to generate the plot for Figure 2A.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "--model-info",
    required=True,
    help="Path to model info csv"
)
parser.add_argument(
    "--bigg-data",
    required=True,
    help="Path to csv with BiGG species"
)
parser.add_argument(
    "--exp-data",
    required=True,
    help="Path to csv with experimental data"
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
from typing import Sequence
from pathlib import Path

import polars as pl
import plotnine as pn

from bates_reconstruction_data import taxonomy
from bates_reconstruction_data.utils import tree


if args.cache is not None:
    TAXDMP = str(Path(args.cache) / "tax.dmp")
    SEGMENTS_DATA = Path(args.cache) / "figure_2a_tree_segments.csv"
    LAYOUT_DATA = Path(args.cache) / "figure_2a_tree_layout.csv"
    DATA_EXISTS = SEGMENTS_DATA.exists() and LAYOUT_DATA.exists()
else:
    TAXDMP = None
    SEGMENTS_DATA = None
    LAYOUT_DATA = None
    DATA_EXISTS = None


def taxonomy_tree(taxa: list[taxonomy.Taxonomy]) -> pl.DataFrame:
    """
    Takes an input list ot taxa and returns a taxonomic tree in the form of a
    DataFrame.

    Returns
    -------
    pl.DataFrame
        A DataFrame with the taxonomic tree. The DataFrame as two columns:
        'node' and 'parent'. Each row represents a node in the tree. The
        'node' column contains the ID of a given node, and the 'parent'
        column contains the ID of that node's parent node.
    """

    tax_levels = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    nodes: dict[str, list] = {"node": [], "parent": []}
    for tx in taxa:
        prev = None
        for level in tax_levels:
            curr = str(getattr(tx, level + "_"))
            if curr not in nodes["node"]:
                nodes["node"].append(curr)
                nodes["parent"].append(prev)
            prev = curr

    return pl.DataFrame(nodes)


def format_species_name(names: Sequence[str]) -> Sequence[str]:
    formatted = []
    for name in names:
        pieces = name.split(maxsplit=2)
        genus = pieces[0]
        sp = pieces[1]

        if sp == "sp.":
            formatted_name = f"$\\it{{{genus}}}$ $\\text{{{sp}}}$"
        else:
            formatted_name = f"$\\it{{{genus} \\ {sp}}}$"

        if len(pieces) > 2:
            formatted_name += f" $\\text{{{pieces[2]}}}$"

        formatted.append(formatted_name)
    return formatted


if args.cache is not None and DATA_EXISTS:
    segments = pl.read_csv(SEGMENTS_DATA)   # ty:ignore[invalid-argument-type]
    layout = pl.read_csv(LAYOUT_DATA)       # ty:ignore[invalid-argument-type]

else:
    # Get a list of bacterial lineages in BiGG
    lineages = taxonomy.get_taxa(
        pl.read_csv(args.bigg_data).get_column("organism"),
        TAXDMP
    )

    # Load selected species
    model_info = pl.read_csv(args.model_info)
    selected = model_info.get_column("species")

    # Get species with experimental data
    species_with_data = (
        pl.read_csv(args.exp_data)
        .join(model_info, on="organism_id")
        .get_column("species")
    )

    # Build the tree
    tree_df = taxonomy_tree(lineages)
    t: dict[str, tree.TreeNode] = tree.build_tree(tree_df, "node", "parent")
    leaf_ids = [n for n in t if len(t[n].children) == 0]

    # Generate a layout of the taxonomic tree
    method = tree.TreeLayout(alpha=1e4, shuffle=False)
    layout = method.layout_nodes({x: t[x] for x in sorted(t)})
    segments = method.node_connections(t, layout)

    # Add metadata
    layout = layout.with_columns(
        selected=pl.col.node.is_in(selected.implode()),
        has_data=pl.col.node.is_in(species_with_data.implode()),
        is_leaf=pl.col.node.is_in(leaf_ids)
    )

    # Save computed layout
    if args.cache is not None:
        segments.write_csv(SEGMENTS_DATA)
        layout.write_csv(LAYOUT_DATA)


# Pull out leaf nodes
leaves = layout.filter(pl.col.is_leaf)

# Plot the tree
plot = (
    pn.ggplot(
        mapping=pn.aes(x="depth", y="width")
    )
    + pn.geom_text(
        pn.aes(label="format_species_name(node)"),
        leaves,
        size=9,
        ha="left",
        nudge_x=2.5,
        show_legend=False,
    )
    + pn.geom_segment(
        pn.aes(xend="depth_end", yend="width_end"),
        segments,
        size=1,
        lineend="projecting"
    )
    + pn.geom_point(
        pn.aes(fill="selected"),
        data=leaves,
        size=2,
        stroke=1,
        color="black",
    )
    + pn.geom_point(
        pn.aes(x="depth+1.5", color="has_data"),
        data=leaves.filter(pl.col.has_data),
        shape="*",
        size=3
    )
    + pn.scale_fill_manual(
        values=["black", "white"],
        limits=[True, False],                   # ty:ignore[invalid-argument-type]
        labels=["Selected", "Not selected"]
    )
    + pn.scale_color_manual(
        values=["#d10202"],
        limits=[True],                          # ty:ignore[invalid-argument-type]
        labels=["Substrate\nutilization data"],
        guide=pn.guide_legend(override_aes={"size": 4})     # ty:ignore[invalid-argument-type]
    )
    + pn.labs(fill="", color="")
    + pn.coord_cartesian(xlim=(0, 32)) # prevent clipping text
    + pn.theme_void()
    + pn.theme(
        figure_size=(2.7, 5),
        legend_position="bottom",
        legend_direction="vertical",
        legend_box_spacing=-0.05,
        legend_text=pn.element_text(size=10)
    )
)

plot.save(args.out)
