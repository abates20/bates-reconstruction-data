from typing import Any, Optional, TypeVar, Iterable
import enum
import random

import polars as pl
import numpy as np


_KeyType = TypeVar("_KeyType")


class ConnectionType(enum.Enum):
    RECTANGULAR = "rectangular"
    SLANTED = "slanted"


class TreeNode:

    def __init__(
        self,
        value: Any,
        parent: Optional["TreeNode"] = None,
        children: Optional[dict[Any, "TreeNode"]] = None,
        metadata: Optional[dict] = None
    ):
        self.value = value
        self.parent = parent
        self.children: dict[Any, "TreeNode"] = children or {}
        self.metadata: dict = metadata or {}

    def __str__(self):
        string = f"TreeNode({repr(self.value)})"

        parent = f"TreeNode({self.parent.value})" if self.parent else None
        string += f"\n\tParent: {parent}"

        elbow = "\u2514"
        hline = "\u2500"
        vline = "\u251c"
        string += "\n\tChildren:"
        for i, child in enumerate(self.children):
            if i == (len(self.children) - 1): # last child
                string += f"\n\t  {elbow}{hline} TreeNode({child})"
            else:
                string += f"\n\t  {vline}{hline} TreeNode({child})"

        string += f"\n\tMetadata: {repr(self.metadata)}"
        return string
    
    def add_child(self, child: "TreeNode"):
        self.children[child.value] = child
    
    def set_parent(self, parent: Optional["TreeNode"]):
        self.parent = parent

    @property
    def descendants(self) -> dict[Any, "TreeNode"]:
        descendants = {}
        descendants.update(self.children)
        for child in self.children.values():
            descendants.update(child.descendants)
        return descendants


def build_tree(data: pl.DataFrame, node_col, parent_col, metadata_cols = None) -> dict[Any, TreeNode]:
    nodes: dict[Any, TreeNode] = {}

    if metadata_cols is None:
        metadata_cols = []

    rows_iter = data.select(node_col, parent_col, *metadata_cols).iter_rows()

    for node_val, parent_val, *metadata in rows_iter:
        default = TreeNode(parent_val) if parent_val is not None else None
        parent = nodes.get(parent_val, default)

        node = nodes.get(node_val, TreeNode(node_val))
        node.metadata = dict(zip(metadata_cols, metadata))
        if parent:
            parent.add_child(node)
            node.set_parent(parent)
            nodes[parent_val] = parent

        nodes[node_val] = node
    
    return nodes


class TreeLayout:

    def __init__(
            self,
            alpha: float = 10.0,
            align_leaves: bool = False,
            shuffle: bool = False,
            random_seed = None
        ):
        self.alpha = alpha
        self.align_leaves = align_leaves
        self.shuffle = shuffle
        self._rng = random.Random(random_seed)

    def layout_nodes(self, tree: dict[Any, TreeNode]) -> pl.DataFrame:
        
        # Get leaf nodes for use later
        leaves = [x for x in tree.values() if len(x.children) == 0]

        # Get depth of each node
        depths = TreeLayout.compute_node_depths(tree)
        if self.align_leaves:
            deepest_leaf = max(leaves, key=lambda x: depths[x.value])
            leaf_depth = depths[deepest_leaf]
            for leaf in leaves:
                depths[leaf] = leaf_depth

        # Width positions
        widths = {x: 0.0 for x in depths}

        # Set widths of leaf nodes (base node is at 0)
        top_nodes = self._sort_nodes(
            n for n in tree.values() if n.parent is None
        )
        self._place_leaves(
            widths,
            top_nodes, 
            list(range(len(top_nodes))), 
            [1.0] * len(top_nodes)
        )

        # Set widths of parent nodes based on position of child nodes
        self._place_rec(widths, leaves)

        # Center widths so that root node is at 0
        root_pos = widths[top_nodes[0].value]
        for node in widths:
            widths[node] = widths[node] - root_pos

        # Return depth and width in dataframe
        node_keys = [k for k in depths]
        layout = (
            pl.DataFrame(
                {
                    "node": node_keys,
                    "depth": [depths[k] for k in node_keys],
                    "width": [widths[k] for k in node_keys]
                }
            )
        )
        return layout
    
    def node_connections(
            self,
            tree: dict[Any, TreeNode],
            layout_data: pl.DataFrame,
            how: ConnectionType = ConnectionType.RECTANGULAR
        ) -> pl.DataFrame:
        how = ConnectionType(how)

        # Dataframe with just position info
        positions = layout_data.select("node", "depth", "width")

        # Get parent for each node
        def get_parent(x):
            node = tree.get(x)
            if node is None:
                return None
            parent = node.parent
            return parent.value if parent is not None else None

        layout_data = (
            layout_data
            .clone()
            .with_columns(
                parent=pl.col.node.map_elements(get_parent, return_dtype=pl.String)
            )
        )

        # Connect parent nodes to children nodes with dataframe join
        segments = (
            layout_data
            .join(
                positions,
                left_on="parent",
                right_on="node",
                suffix="_end"
            )
        )

        if how == ConnectionType.RECTANGULAR:
            segments = segments.with_columns(width_end="width")

            # Add crossing segments
            cross_segs = (
                layout_data
                .group_by("parent")
                .agg(
                    width=pl.col.width.min(),
                    width_end=pl.col.width.max()
                )
                .join(
                    positions.drop("width"),
                    left_on="parent",
                    right_on="node"
                )
                .with_columns(depth_end="depth")
            )

            segments = pl.concat([segments, cross_segs], how="diagonal")
        
        return segments

    def _place_rec(self, widths: dict[Any, float], children: list[TreeNode]):
        parent_dict = {} # Build dict to make sure you get each parent only once
        for child in children:
            if child.parent is not None:
                parent_dict[child.parent.value] = child.parent
        parents = list(parent_dict.values())

        if len(parents) == 0:
            return widths

        for node in parents:
            child_pos = [widths[x] for x in node.children]
            
            # Skip node if any of its children have not been placed yet
            if any(x is None for x in child_pos):
                continue

            widths[node.value] = (min(child_pos) +  max(child_pos)) / 2
        
        return self._place_rec(widths, parents)

    def _place_leaves(
            self,
            width_data: dict[Any, float],
            nodes: list[TreeNode],
            positions: list[float],
            widths: list[float]
        ):
        for node, pos, width in zip(nodes, positions, widths):
            children = self._sort_nodes(node.children.values())

            if len(children) == 0:
                width_data[node.value] = pos
                continue
            
            if len(children) == 1:
                child_positions = [pos]
                child_widths = [width]
            else:
                breaks = self._get_breaks(children, pos, width)
                child_positions = (breaks + np.roll(breaks, -1)) / 2
                child_positions = child_positions.tolist()[:-1] # Take off the extra value at the end
                child_widths = (np.roll(breaks, -1) - breaks).tolist()[:-1]

            self._place_leaves(width_data, children, child_positions, child_widths)

    def _get_breaks(self, nodes: list[TreeNode], parent_y: float, width: float):
        weights = np.ones(len(nodes))
        for i, node in enumerate(nodes):
            n_leaves = len([x for x in node.descendants.values() if len(x.children) == 0])
            weights[i] = weights[i] + (max(n_leaves, 1) * self.alpha)
        
        # normalize weights
        weights /= sum(weights)

        # cumsum to get relative breaks (and stick 0 at beginning)
        relative_breaks = np.insert(weights.cumsum(), 0, 0)

        # Adjust breaks to absolute position
        return width * (relative_breaks - 0.5) + parent_y

    def _sort_nodes(self, nodes: Iterable[TreeNode]) -> list[TreeNode]:
        if self.shuffle:
            nodes = list(nodes)
            self._rng.shuffle(nodes)
            return nodes
        return sorted(nodes, key=lambda x: -len(x.descendants))

    @staticmethod
    def compute_node_depths(tree: dict[_KeyType, TreeNode]) -> dict[_KeyType, float]:
        depths = {}
        for node in tree.values():
            if node.parent is None:
                TreeLayout._compute_depths_rec(depths, node, 0)
        return depths
    
    @staticmethod
    def _compute_depths_rec(depths: dict[Any, float], node: TreeNode, curr_depth: float):
        depths[node.value] = curr_depth
        for child in node.children.values():
            TreeLayout._compute_depths_rec(depths, child, curr_depth+1)
