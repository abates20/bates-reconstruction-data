"""
get_unique_annotations.py
-------------------------

Takes a gene -> rxn map json file and creates a new json file that only has
unique annotations.

In this case, an "annotation" refers to a set of reactions that a gene is mapped
to.
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
args = parser.parse_args()


# Additional imports
import json
import gzip


def open_map(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


# Load input gene -> rxn map
print(f"Loading input gene -> rxn map: {args.input_map}...")
with open_map(args.input_map, "r") as f:
    gene_rxn_map: dict[str, list[str]] = json.load(f)

print(f"Input gene annotation map has {len(gene_rxn_map)} annotations.")

# Identify the unique annotations
unique_annotations = set()
new_map = {}
for gene in sorted(gene_rxn_map):
    annotation = frozenset(gene_rxn_map[gene])
    
    # Skip annotation if already seen
    if annotation in unique_annotations:
        continue
    
    # Otherwise add it to the list of seen ones and to the new map
    unique_annotations.add(annotation)
    new_map[gene] = list(annotation)

print(f"Found {len(new_map)} unique annotations.")

# Save the new set of annotations
with open_map(args.out, "wt") as out:
    json.dump(new_map, out, indent=4)
print(f"Saved unique annotations to {args.out}.")
