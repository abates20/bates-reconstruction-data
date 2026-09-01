"""
calculate_set_metrics.py
------------------------

Calculate metrics (like coverage and diversity) for each resource set based on
the species that were selected.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("input", type=Path, help="Path to resource set data csv")
parser.add_argument("--out", type=Path, required=True, help="Path to output csv with resource set metrics")
parser.add_argument("--reconstruction-dir", type=Path, required=True, help="Path to reconstruction directory")
parser.add_argument("--bigg-info", type=Path, required=True, help="Path to csv with BiGG model info")
args = parser.parse_args()


# Additional imports
import json

import polars as pl

from bates_reconstruction_data import taxonomy
from bates_reconstruction_data.utils import progress


# Load the input resource data
resource_sets = list(
    pl.read_csv(
        args.input,
        schema_overrides={"source_id": pl.String}
    )
    .select("source_id", "fraction")
    .iter_rows()
)

# Load the reference BiGG species and map between model IDs and organism names
bigg_info = pl.read_csv(args.bigg_info)
bigg_species = {
    id_: org for id_, org in bigg_info.select("model_id", "organism").iter_rows()
}
ref_species = list(bigg_species.values())

# Dictionary to hold metrics
metrics_dict = {
    "source_id": [],
    "fraction": [],
    "bigg_coverage": [],
    "diversity": [],
    "selected_species": []
}

# Iterate through resource sets and calculate metrics
for set_id, fraction in progress.track(resource_sets):

    metrics_dict["source_id"].append(set_id)
    metrics_dict["fraction"].append(fraction)

    # Load the selected species
    selected_species_path: Path = args.reconstruction_dir / set_id / "selected_species.txt"
    with selected_species_path.open("r") as f:
        model_ids = f.read().strip().split("\n")
    species = [bigg_species[id_] for id_ in model_ids]
    
    # Calculate metrics
    coverage = taxonomy.weighted_coverage(species, ref_species)
    metrics_dict["bigg_coverage"].append(coverage)

    diversity = taxonomy.taxonomic_diversity(species)
    metrics_dict["diversity"].append(diversity)

    # Add list of selected species for reference
    clean = lambda x: " ".join(str.split(x, maxsplit=2)[:2])
    metrics_dict["selected_species"].append(",".join(clean(x) for x in species))

# Save metrics
pl.DataFrame(metrics_dict).write_csv(args.out)
