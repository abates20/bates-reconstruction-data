"""
prep_organism_info.py
---------------------

Determine gram stain and random seed for each organism.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument(
    "input",
    help="Path to input csv with model info"
)
parser.add_argument(
    "--extra",
    help="Path to a csv with extra data to add"
)
parser.add_argument(
    "--out",
    required=True,
    help="Path to output csv"
)
args = parser.parse_args()


# Additional imports
import polars as pl

from bates_reconstruction_data.taxonomy import gram_stain


# Load model info
model_info = pl.read_csv(args.input)

# Add gram stains
stains = []
for organism_id, species in model_info.select("organism_id", "species").iter_rows():
    stain = gram_stain(species)
    if stain is not None:
        stain = stain.value
    stains.append({"organism_id": organism_id, "gram_stain": stain})

info = (
    model_info
    .select("organism_id", "refseq_accession", "bigg_model_id")
    .join(pl.DataFrame(stains), on="organism_id")
)

# Add extra data if any
if args.extra is not None:
    info = info.join(pl.read_csv(args.extra), on="organism_id")

# Save info
info.write_csv(args.out)
