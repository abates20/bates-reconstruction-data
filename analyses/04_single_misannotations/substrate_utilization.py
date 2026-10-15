"""
substrate_utilization.py
------------------------

This script simulates growth on provided substrates for the given models.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("organism_id", help="ID of bacteria to run predictions for")
parser.add_argument("--out", required=True, type=Path, help="Path to save output data")
parser.add_argument("--set-info", required=True, type=Path, help="Path to csv with annotation set info")
parser.add_argument("--models", required=True, type=Path, help="Path to models directory")
parser.add_argument("--exp-data", required=True, type=Path, help="Path to csv with experimental growth data")
parser.add_argument("--mediadb", required=True, type=Path, help="Path to medium database (as .tsv file)")
parser.add_argument("--processes", type=int, default=None, help="Number of parallel processes")
parser.add_argument("--solver", default=None, help="What solver to use")
parser.add_argument("--log", default="INFO", help="Logging level")
args = parser.parse_args()


# Additional imports
from itertools import repeat, chain
import logging

import polars as pl
import cobra

from bates_reconstruction_data.utils import parallel


# Global params
MEDIA = "minimal_noCarbon[bigg]"
SUBSTRATE_UPTAKE_BOUND = 10.0
GROWTH_CUTOFF = 1e-6

# Configure logging
logging.basicConfig(
    format="[%(asctime)s] : %(levelname)s : %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=args.log
)

# Set cobra solver
COBRA_CONFIG = cobra.Configuration()
if args.solver is not None:
    COBRA_CONFIG.solver = args.solver
logging.info(f"Cobra solver set to {COBRA_CONFIG.solver}")


def get_medium(name: str, db: pl.DataFrame) -> dict[str, float]:
    """
    Get a medium specification (as a dict) from the media db.
    """
    medium = {}
    for cpd in db.filter(pl.col.medium == name).get_column("compound"):
        medium[f"EX_{cpd}_e"] = 1000.0
    return medium


def set_medium(
    model: cobra.Model,
    medium: dict[str, float]
):
    """
    Set the medium for a model from the medium dict.

    Ignores exchange reactions if they do not exist in the model.
    """
    clean_medium = {
        ex: medium[ex] for ex in medium if model.reactions.has_id(ex)
    }
    model.medium = clean_medium


def run_predictions(
    model_path: Path,
    set_id: str,
    substrates: list[tuple[str, list[str]]],
    base_medium: dict[str, float],
) -> list[dict]:
    """
    Run substrate utilization predictions for the specified model.
    """
    TOTAL = len(substrates)
    logging.info(f"Running {TOTAL} predictions for set_id={set_id}")

    model = cobra.io.load_json_model(model_path)
    results = []

    for i, (compound, bigg_ids) in enumerate(substrates, start=1):
        test_medium = base_medium.copy()

        in_model = False
        for bigg_id in bigg_ids:
            if model.metabolites.has_id(f"{bigg_id}_e"):
                in_model = True
                test_medium[f"EX_{bigg_id}_e"] = SUBSTRATE_UPTAKE_BOUND

        if in_model:
            with model as m:
                set_medium(m, test_medium)
                growth = m.slim_optimize()
        else:
            growth = 0.0

        results.append(
            {
                "set_id": set_id,
                "compound": compound,
                "growth": growth > GROWTH_CUTOFF
            }
        )

        if i % 10 == 0 or i == TOTAL:
            logging.info(f"[set_id={set_id}] Prediction progress: {i} / {TOTAL} ({i / TOTAL:.1%})")

    return results


# Load the info about annotation sets
annotation_sets = pl.read_csv(
    args.set_info,
    schema_overrides={"set_id": pl.String}
)
logging.debug(f"Annotation sets:\n{annotation_sets}")

# Load the media DB and get the compounds in the media
medium_df = pl.read_csv(args.mediadb, separator="\t")
medium = get_medium(MEDIA, medium_df)

# Load the growth data
growth_data = pl.read_csv(args.exp_data)

# Get substrates from growth data
substrates = list(
    growth_data
    .drop_nulls("bigg_id")
    .unique("compound")
    .select("compound", pl.col.bigg_id.str.split(";"))
    .iter_rows()
)

# Run predictions for the reconstructed models
set_ids = list(annotation_sets.get_column("set_id"))
results = parallel.apply(
    run_predictions,
    [args.models / f"{x}.json" for x in set_ids],
    set_ids,
    repeat(substrates),
    repeat(medium),
    max_workers=args.processes,
    executor_type="process",
    show_progress=False
)

# Save results
(
    pl.DataFrame(chain.from_iterable(results))
    .select(
        pl.lit(args.organism_id).alias("organism_id"),
        pl.all()
    )
    .write_csv(args.out)
)