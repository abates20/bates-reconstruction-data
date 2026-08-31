"""
reaction_matrix.py
------------------

Script to go through a set of reconstructed models and create a reaction
matrix of 0s and 1s to record which reactions are present/absent in each model.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("--model-info", required=True, type=Path, help="Path to csv with model info")
parser.add_argument("--set-info", required=True, type=Path, help="Path to file with resource set info")
parser.add_argument("--models", required=True, type=Path, help="Path to base models directory")
parser.add_argument("--out", type=Path, required=True, help="Path to output reaction matrix")
parser.add_argument("--include-bigg", action="store_true", help="Include rows for curated BiGG models")
parser.add_argument("--to-mnx", action="store_true", help="Map reaction IDs to the MetaNetX namespace")
parser.add_argument("--xref-cache", type=Path, default=None, help="Path to cache reaction xref table (optional)")
args = parser.parse_args()


# Additional imports
from typing import Optional
import shutil
import json
import io

import requests
import polars as pl
import cobra

from bates_reconstruction_data.utils import progress
from bates_reconstruction_data import bigg


# Load the reac_xref table (for mapping reactions between ModelSEED and BiGG) from MetaNetX
if args.to_mnx:
    reac_xref_path: Optional[Path] = args.xref_cache

    if reac_xref_path is not None and reac_xref_path.exists():
        reac_xref = pl.read_csv(reac_xref_path)

    else:
        # Download and clean it up if it hasn't already been downloaded
        print("Downloading and cleaning up MetaNetX xref table to map reaction IDs...")
        
        reac_xref_url = "https://www.metanetx.org/ftp/4.5/reac_xref.tsv"
        response: requests.Response = requests.get(reac_xref_url)
        response.raise_for_status()

        reac_xref = pl.read_csv(
            io.BytesIO(response.content),
            separator="\t",
            comment_prefix="#",
            quote_char=None,
            has_header=False,
            new_columns=["xref", "mnx_id", "description"]
        )
        
        # Filter to just BiGG and ModelSEED reactions with some reformatting of the data
        reac_xref = (
            reac_xref
            .with_columns(
                pl.col.xref
                .str.split_exact(":", 1)
                .struct.rename_fields(["xref_prefix", "xref_id"])
            )
            .unnest("xref")
            .filter(
                pl.col.xref_prefix.str.contains("(bigg|seed).reaction"),
                ~pl.col.xref_id.str.contains("^(R|EX|SK|DM)_"),
                pl.col.mnx_id.ne("EMPTY")
            )
            .with_columns(
                pl.col.xref_prefix
                .str.split(".")
                .list.first().alias("xref_source")
            )
            .select("xref_source", "xref_id", "mnx_id")
        )

        # Just include metanetx IDs that map to both BiGG and ModelSEED reactions
        mnx_ids_to_keep = (
            reac_xref
            .group_by("mnx_id")
            .agg(pl.col.xref_source.unique())
            .filter(pl.col.xref_source.list.len() > 1)
            .select("mnx_id")
        )
        reac_xref = reac_xref.join(mnx_ids_to_keep, on="mnx_id", how="semi")

        # Save cleaned up table
        if reac_xref_path is not None:
            reac_xref.write_csv(reac_xref_path)
            print(f"Saved xref table to {reac_xref_path}.")
    

# Load the model info
model_info = pl.read_csv(args.model_info)
models = list(
    model_info
    .select("organism_id", "refseq_accession")
    .iter_rows()
)


# Load the info about reconstruction resource sets
resource_sets = list(
    pl.read_csv(
        args.set_info,
        schema_overrides={"source_id": pl.String}
    )
    .get_column("source_id")
)


# Compile the data about which reactions are present/absent in each model
rxn_data_pieces: list[pl.DataFrame] = []

for source_id in resource_sets:
    print(f"Compiling reaction data for source_id={repr(source_id)} models:")

    for org_id, accession in progress.track(models):
        model = cobra.io.read_sbml_model(args.models / source_id / f"{accession}.sbml")
        df = pl.DataFrame(
            {
                "organism_id": org_id,
                "source_id": source_id,
                "rxn_id": [r.id for r in model.reactions],
                "is_present": 1
            }
        )
        rxn_data_pieces.append(df)


# Add info for BiGG models if specified
if args.include_bigg:
    print("Compiling reaction data for curated BiGG models...")

    for org_id, model_id in model_info.select("organism_id", "bigg_model_id").iter_rows():
        model = bigg.get_model(model_id)
        df = pl.DataFrame(
            {
                "organism_id": org_id,
                "source_id": "BiGG",
                "rxn_id": [r["id"] for r in model["reactions"]],
                "is_present": 1
            }
        )
        rxn_data_pieces.append(df)


# Combine all the data into one dataframe
rxn_data = pl.concat(rxn_data_pieces, how="diagonal")

# Map to MetaNetX IDs (if specified)
if args.to_mnx:
    rxn_data = rxn_data.join(
        reac_xref,
        left_on="rxn_id",
        right_on="xref_id",
        how="inner",
        validate="m:1"
    )
    pivot_col = "mnx_id"
else:
    pivot_col = "rxn_id"

# Clean up and reformat data
rxn_data = (
    rxn_data
    .unique(subset=["organism_id", "source_id", pivot_col])
    .pivot(
        pivot_col,
        values="is_present",
        index=["organism_id", "source_id"]
    )
    .fill_null(0)
)
        
# Save the data
rxn_data.write_csv(args.out)
print(f"Saved the reaction data to {args.out}.")
