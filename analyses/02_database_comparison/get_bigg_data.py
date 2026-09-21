"""
get_bigg_data.py
----------------

Download data from BiGG about organisms, reactions, and genes to compare against
other databases.
"""


from argparse import ArgumentParser
from pathlib import Path

parser = ArgumentParser()
parser.add_argument("--counts-out", required=True, help="Path to output counts")
parser.add_argument("--functions-out", required=True, help="Path to output list of functions")
parser.add_argument("-t", "--threads", type=int, default=5, help="Number of threads for parallel downloading")
parser.add_argument("--log", default="INFO", help="Logging level")
args = parser.parse_args()


# Additional imports
import itertools
import re
import logging

import polars as pl
import requests

from bates_reconstruction_data.taxonomy import Taxonomy, NoTaxonomyDataError
from bates_reconstruction_data import bigg
from bates_reconstruction_data.utils import parallel


# Configure logging
logging.basicConfig(
    format="[%(asctime)s] : %(levelname)s : %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=args.log
)


def split_gpr(gpr: str) -> set[str]:
    not_genes = r"\s?(AND|OR|[\(\)])\s?"
    gpr = re.sub(not_genes, " ", gpr, flags=re.IGNORECASE)
    return set(gpr.split())


# Download the info about all the BiGG models
logging.info("Getting model metadata from BiGG...")
response = bigg.get_model_metadata()
model_info = pl.from_dicts(response.get("results"))

# Add some taxonomy info
logging.info("Adding taxononmy information from NCBI...")
model_info = (
    model_info
    .with_columns(
        species=(
            pl.col.organism
            .str.split(" ")
            .list.slice(0, 2)
            .list.join(" ")
        ) # Add species column
    )
)

tax_data = {
    "species": [],
    "domain": []
}
for species in model_info.get_column("species").unique():
    try:
        tax = Taxonomy.from_name(species)
    except NoTaxonomyDataError:
        print(f"No taxonomy info for {species}. Skipping...")
    else:
        tax_data["species"].append(species)
        tax_data["domain"].append(tax.domain_)

model_info = (
    model_info
    .join(
        pl.DataFrame(tax_data),
        on="species"
    )
    .filter(pl.col.domain == "Bacteria") # just keep bacteria species
)

# Download models to get genes and reactions
logging.info("Downloading and processing bacteria models from BiGG...")
with requests.Session() as s:
    models = parallel.apply(
        bigg.get_model,
        model_info.get_column("bigg_id"),
        itertools.repeat(s),
        threads=args.threads
    )

# compile reactions and associated genes from models
gene_rxn_dct = {
    "bigg_id": [],
    "gene": [],
    "reaction": [],
    "ec_number": []
}
for model in models:
    for rxn in model["reactions"]:
        gpr = rxn["gene_reaction_rule"]
        genes = [None] if gpr == "" else split_gpr(gpr)

        for gene in genes:
            gene_rxn_dct["bigg_id"].append(model["id"])
            gene_rxn_dct["gene"].append(gene)
            gene_rxn_dct["reaction"].append(rxn["id"])
            gene_rxn_dct["ec_number"].append(rxn["annotation"].get("ec-code"))

model_info = (
    model_info.join(
        pl.DataFrame(gene_rxn_dct),
        on="bigg_id"
    )
    .select("bigg_id", "species", "gene", "reaction", "ec_number")
    .explode("ec_number")
)
logging.debug(str(model_info))

# Save bacteria-associated functions
(
    model_info
    .filter(pl.col.gene.is_not_null())
    .select(
        gene=pl.col.bigg_id + "." + pl.col.gene,
        reaction="reaction"
    )
    .write_csv(args.functions_out)
)

# Save counts per organsim (i.e., per model)
bigg_counts = (
    model_info
    .filter(pl.col.gene.is_not_null())
    .group_by("bigg_id", "species")
    .agg(
        n_genes=pl.col.gene.unique().len(),
        n_rxns=pl.col.reaction.unique().len()
    )
)
bigg_counts.write_csv(args.counts_out)
