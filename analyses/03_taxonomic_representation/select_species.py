"""
select_species.py
-----------------

Randomly selects a specific fraction of the supplied BiGG models and writes
their IDs to the specified output file.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("input", type=Path, help="Path to csv with BiGG model info")
parser.add_argument("--fraction", type=float, required=True, help="Fraction of models to select")
parser.add_argument("--seed", type=int, help="Random seed")
parser.add_argument("--out", type=Path, help="Output path to save list of selected model IDs")
args = parser.parse_args()


# Additional imports
import polars as pl


# Load input model info
models = pl.read_csv(args.input)

# Select models
selected = models.sample(fraction=args.fraction, seed=args.seed)

# Save list of model IDs
with open(args.out, "w") as f:
    for model_id in selected.get_column("model_id").unique():
        f.write(f"{model_id}\n")
