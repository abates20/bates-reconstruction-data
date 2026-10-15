"""
process_reactions.py
--------------------

Compile a matrix of what reactions are present in a set of models.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("organism_id", help="ID of bacteria to run predictions for")
parser.add_argument("--out", required=True, type=Path, help="Path to save output data")
parser.add_argument("--set-info", required=True, type=Path, help="Path to csv with annotation set info")
parser.add_argument("--models", required=True, type=Path, help="Path to models directory")
parser.add_argument("--include-bigg", help="Add reactions from curated BiGG model")
parser.add_argument("--log", default="INFO", help="Logging level")
args = parser.parse_args()


# Additional imports
import json
import logging

import polars as pl

from bates_reconstruction_data import bigg


# Configure logging
logging.basicConfig(
    format="[%(asctime)s] : %(levelname)s : %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=args.log
)


# Load the info about annotation sets
annotation_sets = pl.read_csv(
    args.set_info,
    schema_overrides={"set_id": pl.String}
)
logging.debug(f"Annotation sets:\n{annotation_sets}")

# Pull out the reaction in each model
logging.info("Compiling reaction data for models")
set_ids = annotation_sets.get_column("set_id")
n_sets = len(set_ids)
reaction_data = []

for i, set_id in enumerate(set_ids, start=1):

    with open(args.models / f"{set_id}.json") as f:
        model = json.load(f)
    
    data = {"set_id": set_id}
    data.update((r["id"], 1) for r in model["reactions"])
    reaction_data.append(data)

    if i % 50 == 0 or i == n_sets:
        logging.info(f"Processing models: {i} / {n_sets} completed")

# Add info for BiGG model if specified
if args.include_bigg is not None:
    logging.info(f"Adding reactions for BiGG model {args.include_bigg}")

    model = bigg.get_model(args.include_bigg)
    data = {"set_id": "BiGG"}
    data.update((r["id"], 1) for r in model["reactions"])
    reaction_data.append(data)

# Save the results
(
    pl.DataFrame(reaction_data, infer_schema_length=None)
    .select(
        pl.lit(args.organism_id).alias("organism_id"),
        pl.all()
    )
    .fill_null(0)
    .write_csv(args.out)
)
