"""
make_gene_rxn_map.py
--------------------

Builds a gene -> rxn map for Reconstructor from a list of BiGG models.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument(
    "models",
    nargs="+",
    help="List of BiGG model IDs"
)
parser.add_argument(
    "--out",
    required=True,
    type=Path,
    help="Path to save output map"
)
parser.add_argument(
    "--threads",
    default=5,
    type=int,
    help="Number of threads for parallel downloading"
)

args = parser.parse_args()


# Additional imports
import json
import re
from collections import defaultdict
import gzip
import itertools

import requests

from bates_reconstruction_data import bigg
from bates_reconstruction_data.utils import parallel


def parse_gpr(gpr: str):
    not_genes = r"\s?(AND|OR|[\(\)])\s?"
    gpr = re.sub(not_genes, " ", gpr, flags=re.IGNORECASE)
    return set(gpr.split())


if __name__ == "__main__":

    # Get the models from BiGG
    with requests.Session() as session:
        models = parallel.apply(
            bigg.get_model,
            args.models,
            itertools.repeat(session),
            max_workers=args.threads
        )

    # Create the gene -> reaction map
    gene_rxn_map: defaultdict[str, list] = defaultdict(list)
    for model_id, model in zip(args.models, models):
        for rxn in model["reactions"]:
            for gene in parse_gpr(rxn["gene_reaction_rule"]):
                gene_key = f"{model_id}.{gene}"
                gene_rxn_map[gene_key].append(rxn["id"])

    # Save the map
    opener = gzip.open if args.out.suffix == ".gz" else open
    with opener(args.out, "wt") as f:
        json.dump(gene_rxn_map, f)
