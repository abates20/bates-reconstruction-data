"""
substrate_utilization.py
------------------------

This script simulates growth on provided substrates for the given models.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("--out", required=True, type=Path, help="Path to save output data")
parser.add_argument("--model-info", required=True, type=Path, help="Path to csv with model info")
parser.add_argument("--set-info", required=True, type=Path, help="Path to file with resource set info")
parser.add_argument("--models", required=True, type=Path, help="Path to base models directory")
parser.add_argument("--data", required=True, type=Path, help="Path to csv with experimental growth data")
parser.add_argument("--media-db", required=True, type=Path, help="Path to medium database (as .tsv file)")
args = parser.parse_args()


# Additional imports
import polars as pl
import cobra

from bates_reconstruction_data.utils import progress


MEDIA = "minimal_noCarbon"
SUBSTRATE_UPTAKE_BOUND = 10.0
GROWTH_CUTOFF = 1e-6


def get_medium(name: str, db: pl.DataFrame):
    medium = {}
    for cpd in db.filter(pl.col.medium == name).get_column("compound"):
        medium[f"EX_{cpd}_e"] = 1000.0
    return medium


def set_medium(
        model: cobra.Model,
        medium: dict[str, float]
    ):
    clean_medium = {
        ex: medium[ex] for ex in medium if model.reactions.has_id(ex)
    }
    model.medium = clean_medium


# Load the model info
model_info = pl.read_csv(args.model_info)

# Load the info about reconstruction resource sets
resource_sets = list(
    pl.read_csv(
        args.set_info,
        schema_overrides={"source_id": pl.String}
    )
    .select("source_id", "namespace")
    .iter_rows()
)

# Load the media DB
medium_df = pl.read_csv(args.media_db, separator="\t")

# Load the growth data
growth_data = pl.read_csv(args.data)

# Prep the growth data
growth_data = (
    growth_data
    .unique(["organism_id", "compound"])
    .join(model_info, on="organism_id")
    .select(
        "organism_id",
        source_id=pl.lit("experimental"),
        compound="compound",
        bigg_id=pl.col.bigg_id.str.split(";"),
        seed_id=pl.col.seed_id.str.split(";"),
        growth="growth",
        in_model=pl.lit(None, dtype=pl.Boolean)
    )
    .drop_nulls(["bigg_id", "seed_id"])
)

# Get substrate list
substrates = (
    growth_data
    .select("compound", "bigg_id", "seed_id")
    .unique()
)

# Dictionary to hold predicted growth data
growth_pred_dict = {
    "organism_id": [],
    "source_id": [],
    "compound": [],
    "bigg_id": [],
    "seed_id": [],
    "growth": [],
    "in_model": []
}

# Run predictions for the reconstructed models
model_info_iter = (
    model_info
    .join(growth_data, on="organism_id", how="semi")
    .select("species", "refseq_accession", "organism_id")
    .iter_rows()
)
for species, accession, org_id in model_info_iter:
    print(f"Running predictions for {species} models (accession: {accession})...")
    
    for source_id, namespace in progress.track(resource_sets):

        # Load model
        model_path = args.models / source_id / f"{accession}.sbml"
        model = cobra.io.read_sbml_model(model_path)
            
        for compound, bigg_ids, seed_ids in substrates.iter_rows():
            test_medium = get_medium(
                MEDIA + ("[bigg]" if namespace == "BiGG" else "[seed]"),
                medium_df
            )
            substrate_ids = bigg_ids if namespace == "BiGG" else seed_ids
            
            in_model = False
            for substrate_id in substrate_ids:
                for compartment in ["c", "e", "p"]:
                    if model.metabolites.has_id(f"{substrate_id}_{compartment}"):
                        in_model = True
                        test_medium[f"EX_{substrate_id}_e"] = SUBSTRATE_UPTAKE_BOUND

            if in_model:
                with model as m:
                    set_medium(m, test_medium)
                    growth = m.slim_optimize()
            else:
                growth = 0.0

            growth_pred_dict["organism_id"].append(org_id)
            growth_pred_dict["source_id"].append(source_id)
            growth_pred_dict["compound"].append(compound)
            growth_pred_dict["bigg_id"].append(bigg_ids)
            growth_pred_dict["seed_id"].append(seed_ids)
            growth_pred_dict["growth"].append(growth > GROWTH_CUTOFF)
            growth_pred_dict["in_model"].append(in_model)


# Save the growth data and predictions together
growth_preds = pl.DataFrame(growth_pred_dict)

growth_data = (
    pl.concat([growth_data, growth_preds])
    .with_columns(
        pl.col.bigg_id.list.join(";"),
        pl.col.seed_id.list.join(";")
    )
)
growth_data.write_csv(args.out)
