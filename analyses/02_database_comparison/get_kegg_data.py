"""
get_kegg_data.py
----------------

Get a data from KEGG including bacteria and the orthologies and the reactions
they are linked to.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("--tax-nodes", required=True, help="Path to NCBI taxonomy nodes dump")
parser.add_argument("--lineages", help="Path to NCBI lineages dump (optional)")
parser.add_argument("--counts-out", required=True, help="Path to output counts")
parser.add_argument("--functions-out", required=True, help="Path to output list of functions")
parser.add_argument("--threads", type=int, default=5, help="Number of threads for parallel downloading")
parser.add_argument("--log", default="INFO", help="Logging level")
args = parser.parse_args()


# Additional imports
import io
import logging

import requests
import polars as pl

from bates_reconstruction_data import taxonomy, kegg
from bates_reconstruction_data.utils import progress


# Configure logging
logging.basicConfig(
    format="[%(asctime)s] : %(levelname)s : %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=args.log
)


def read_dmp(file, columns) -> pl.DataFrame:
    with open(file) as f:
        df = pl.read_csv(
            io.StringIO(f.read().replace("\t|\t", "\0").replace("\t|\n", "\n")),
            separator="\0",
            eol_char="\n",
            quote_char=None,
            has_header=False,
            new_columns=columns
        )
    return df


# Load taxonomy nodes
node_columns = [
    "tax_id",
    "parent_tax_id",
    "rank",
    "embl_code",
    "division_id",
    "inherited_div_flag",
    "genetic_code_id",
    "inherited_GC_flag",
    "mitochondrial_genetic_code_id",
    "inherited_MGC_flag",
    "GenBank_hidden_flag",
    "hidden_subtree_root_flag",
    "comments"
]
nodes = read_dmp(args.tax_nodes, node_columns)

# Load NCBI bacteria lineages
ncbi_lineages = (
    taxonomy.get_taxonomy_data(args.lineages)
    .filter(pl.col.domain == "Bacteria")
    .join(
        nodes.select("tax_id", "rank"),
        on="tax_id",
        how="inner"
    )
    .unique("tax_id")
)

# Get KEGG bacteria with lineages
logging.info("Querying KEGG for organisms...")
response: requests.Response = kegg.request("https://rest.kegg.jp/link/taxonomy/genome")
response.raise_for_status()

kegg_orgs = (
    pl.read_csv(
        io.BytesIO(response.content),
        has_header=False,
        new_columns=["genome_id", "tax_id"],
        separator="\t"
    )
    .select(
        kegg_org=pl.col.genome_id.str.strip_prefix("gn:"),
        tax_id=pl.col.tax_id.str.extract(r"taxid:(\d+)").cast(pl.Int64))
    .unique("kegg_org")
    .join(ncbi_lineages, on="tax_id")
    .with_columns(
        species=(
            pl.when(pl.col.rank == "species")
            .then("tax_name")
            .otherwise("species")
        )
    )
    .select("kegg_org", "species")
)


# Get the list of orthologies with reactions
logging.info("Querying KEGG for KO -> reaction relationships...")
response: requests.Response = kegg.request("https://rest.kegg.jp/link/reaction/ko")
response.raise_for_status()

ko_rxn_map = (
    pl.read_csv(
        io.BytesIO(response.content),
        separator="\t",
        has_header=False,
        new_columns=["ko", "reaction"],
        quote_char=None
    )
    .select(
        pl.col.ko.str.strip_prefix("ko:"),
        pl.col.reaction.str.strip_prefix("rn:")
    )
)
logging.debug(str(ko_rxn_map))


orthologies: list[str] = (
    ko_rxn_map
    .get_column("ko")
    .unique()
    .to_list()
)

# Retreive each orthology
logging.info("Downloading the info for each orthology...")
ko_results = kegg.batched_request(
    "https://rest.kegg.jp/get/{}",
    orthologies,
    threads=args.threads
)

# Prepare KO -> gene map
gene_map = {
    "ko": [],
    "gene": []
}

logging.info("Parsing orthology data...")
for ko_result in progress.track(ko_results):
    try:
        ko_result.raise_for_status()
    except requests.HTTPError:
        logging.warning("Error downloading KO information", exc_info=True)
    else:
        for ko_json in kegg.FlatFileReader.parse(ko_result.text):

            ko_id = ko_json.get("ENTRY", {}).get("entry_id")
            if ko_id is None:
                continue # Skip orthology if it's missing an ID

            genes = [x["org_genes"][0]["gene_id"] for x in ko_json.get("GENES", [])]

            if len(genes) == 0:
                continue # Skip orthology if it doesn't include any genes

            gene_map["ko"].extend([ko_id] * len(genes))
            gene_map["gene"].extend(genes)


#  Filter to just bacteria genes
kegg_genes = (
    pl.DataFrame(gene_map)
    .with_columns(kegg_org=pl.col.gene.str.split(":").list.first())
    .join(kegg_orgs, on="kegg_org")
)
logging.debug(str(kegg_genes))

# Save bacteria-associated functions
(
    ko_rxn_map
    .join(kegg_genes, on="ko", how="semi")
    .write_csv(args.functions_out)
)

# Save counts per organsim
kegg_counts = (
    kegg_genes
    .join(ko_rxn_map, on="ko")
    .group_by("kegg_org", "species")
    .agg(
        n_genes=pl.col.gene.unique().len(),
        n_rxns=pl.col.reaction.unique().len()
    )
)
kegg_counts.write_csv(args.counts_out)
