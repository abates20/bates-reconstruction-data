"""
define_annotation_sets.py
-----------------------

Define all of the annotation sets for a particular organism. There are 'n'
annotation sets per gene from the list of blast hits, and each annotation set
has a random function assigned to the selected gene. Creates a table specifying
which gene is selected in each annotation set and what the new function is for
that gene in that set.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument(
    "blast_results",
    type=Path,
    help="Path to blast results file"
)
parser.add_argument(
    "-o", "--out",
    required=True,
    type=Path,
    help="Path to where output table will be saved"
)
parser.add_argument(
    "--annotation-pool",
    required=True,
    type=Path,
    help="Path to json containing the pool of annotations that can be pulled from."
)
parser.add_argument(
    "-n", "--sets-per-gene",
    type=int,
    default=3,
    help="Number of sets to create for each gene (default=3)"
)
parser.add_argument(
    "--seed",
    type=int,
    default=None,
    help="Optional random seed"
)
args = parser.parse_args()


# Additional imports
import json
import random
import gzip

import polars as pl


def load_blast(path: Path) -> pl.DataFrame:
    """
    Loads a blast output file into a DataFrame.
    """
    default_diamond_cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore"
    ]
    return pl.read_csv(
        path,
        has_header=False,
        new_columns=default_diamond_cols,
        separator="\t"
    )


def open_map(path: Path, mode: str):
    """
    Opens a json file (that may be gzipped).
    """
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


# Load input blast hits
print(f"Loading input blast hits: {args.blast_results}...")
blast_results = load_blast(args.blast_results)

# Load the annotation pool
print(f"Loading annotation pool: {args.annotation_pool}...")
with open_map(args.annotation_pool, "r") as f:
    annotation_pool = json.load(f)
genes_pool = list(sorted(annotation_pool))

# Set the random seed
if args.seed is not None:
    random.seed(args.seed)

# Define annotation sets with randomized functions
print(f"Defining annotation sets...")
annotation_sets = [{"set_id": "0000", "gene": None, "new_function": None}]
genes = blast_results.get_column("sseqid").unique().sort()
i = 1
for gene in genes:
    for _ in range(args.sets_per_gene):
        annotation_sets.append(
            {
                "set_id": f"{i:04d}",
                "gene": gene,
                "new_function": annotation_pool[random.choice(genes_pool)]
            }
        )
        i += 1

annotation_set_df = (
    pl.DataFrame(annotation_sets)
    .with_columns(pl.col.new_function.list.join(";"))
)

# Save annotation set specifications
annotation_set_df.write_csv(args.out)
print(f"Saved annotation sets info to {args.out}.")
