"""
randomize_annotations.py
------------------------

Randomly mix up a fraction of annotations from a set of gene annotations.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument(
    "input_map",
    type=Path,
    help="Path to input gene -> reaction map"
)
parser.add_argument(
    "-o", "--out",
    required=True, 
    type=Path,
    help="Path to where new map will be saved"
)
parser.add_argument(
    "--annotation-pool",
    required=True,
    type=Path,
    help="Path to json containing the pool of annotations that can be pulled from."
)
parser.add_argument(
    "--fraction",
    type=float,
    help="Fraction of annotations to randomize"
)
parser.add_argument(
    "--seed",
    type=int,
    default=None,
    help="Optional random seed"
)
parser.add_argument(
    "--log",
    type=Path,
    default=None,
    help="Path to log file"
)
args = parser.parse_args()


# Additional imports
from typing import Optional
import json
import random
import gzip
import os
from contextlib import redirect_stdout

import polars as pl


LOG= Path(args.log) if args.log is not None else os.devnull


def randomize_annotations(gene_rxn_map: dict[str, list[str]], pool: dict[str, list[str]], fraction: float, seed: Optional[int]) -> dict[str, list[str]]:
    # Set the random seed
    if seed is not None:
        random.seed(seed)
    
    new_map = gene_rxn_map.copy()

    # Select genes to change their annotations
    k = round(fraction * len(gene_rxn_map))
    genes = random.sample(sorted(gene_rxn_map), k)
    
    # Randomly assign different annotations
    new_labels = random.choices(sorted(pool), k=k)
    for gene, new_label in zip(genes, new_labels):
        print(f"\t{gene:30}: {gene_rxn_map[gene]} -> {pool[new_label]}")
        new_map[gene] = pool[new_label]
    
    return new_map


def open_map(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


# Load input gene -> rxn map
print(f"Loading input gene -> rxn map: {args.input_map}...")
with open_map(args.input_map, "r") as f:
    gene_rxn_map = json.load(f)

# Load the annotation pool
print(f"Loading annotation pool: {args.annotation_pool}...")
with open_map(args.annotation_pool, "r") as f:
    annotation_pool = json.load(f)

# Randomize annotations
print("Generating randomized annotations...")
with open(LOG, "a") as f, redirect_stdout(f):
    print(f"Randomly assigning new annotations to {args.fraction:.1%} selected genes:")
    new_map = randomize_annotations(gene_rxn_map, annotation_pool, args.fraction, args.seed)

# Save new map
with open_map(args.out, "wt") as out:
    json.dump(new_map, out)
print(f"Saved new gene -> rxn map to {args.out}.")
