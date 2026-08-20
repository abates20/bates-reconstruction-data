"""
prep_GPR_files.py
-----------------

Prepares files containing the relationships between genes and reactions for
Reconstructor and CarveMe.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("--bigg-gprs", type=Path, required=True, help="Path to output BiGG GPR table for CarveMe")
parser.add_argument("--bigg-map", type=Path, required=True, help="Path to output BiGG gene-rxn map for Reconstructor")
parser.add_argument("--kegg-gprs", type=Path, required=True, help="Path to output KEGG GPR table for CarveMe")
parser.add_argument("--kegg-map", type=Path, required=True, help="Path to output KEGG gene-rxn map for Reconstructor")
args = parser.parse_args()


# Additional imports
import gzip
import json

import polars as pl

from bates_reconstruction_data.reconstructor import resources
from bates_reconstruction_data.utils import download


# Download BiGG GPR file from CarveMe for creating a gene -> reaction map
url = "https://github.com/cdanielmachado/carveme/raw/refs/heads/master/carveme/data/generated/bigg_gprs.csv.gz"
download.download(url, args.bigg_gprs)

# Create the BiGG gene -> reaction map
with gzip.open(args.bigg_gprs) as f:
    bigg_gprs = pl.read_csv(f) # ty:ignore[invalid-argument-type]

bigg_gprs = (
    bigg_gprs
    .select(
        pl.col("gene").str.strip_prefix("G_"),
        pl.col("reaction").str.strip_prefix("R_"),
        "model"
    ).with_columns(
        bigg_gene=pl.col("model") + "." + pl.col("gene")
    )
)

bigg_gene_rxn_map: dict[str, list[str]] = {}
for gene, rxns in bigg_gprs.group_by("bigg_gene").agg("reaction").iter_rows():
    bigg_gene_rxn_map[gene] = rxns

with gzip.open(args.bigg_map, "wt") as f:
    json.dump(bigg_gene_rxn_map, f)


# Load the Reconstructor gene -> rxn map
kegg_gene_rxn_map = resources.get_gene_rxn_map()

# Save copy of map to output location (with gene IDs updated to match DIAMOND db)
with gzip.open(args.kegg_map, "wt") as f:
    json.dump(
        {k.replace(":", "."): v for k, v in kegg_gene_rxn_map.items()},
        f
    )

# Create GPR table
gpr_dict = {
    "gene": [],
    "protein": [],
    "reaction": [],
    "model": []
}
for gene, rxns in kegg_gene_rxn_map.items():
    organism, gene_id = gene.split(":")
    for rxn_id in rxns:
        gpr_dict["gene"].append(f"G_{gene_id}")
        gpr_dict["protein"].append(f"P_{gene_id}")
        gpr_dict["reaction"].append(f"R_{rxn_id}")
        gpr_dict["model"].append(organism)
gprs = pl.DataFrame(gpr_dict)

# Save GPR data
if str(args.kegg_gprs).endswith(".gz"):
    opener = gzip.open
else:
    opener = open
with opener(args.kegg_gprs, "wt") as f:
    gprs.write_csv(f)
