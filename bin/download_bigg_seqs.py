"""
download_bigg_seqs.py
---------------------

Downloads the amino acid sequences for all genes in a BiGG model and saves the
sequences to a fasta file.
"""


from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument("model_id", help="BiGGG model ID")
parser.add_argument("--out", help="Path to output fasta file")
parser.add_argument("--refseq-accession", help="RefSeq accession to retreive missing sequences from NCBI")
parser.add_argument("--threads", type=int, default=5, help="Number of threads for parallel downloading")
args = parser.parse_args()


# Additional imports
from typing import Optional
import itertools
from tempfile import TemporaryDirectory
import os

import requests
import polars as pl

from bates_reconstruction_data import bigg
from bates_reconstruction_data.ncbi_genomes import NCBI, AnnotationType
from bates_reconstruction_data.utils import parallel


def read_gff(file) -> pl.DataFrame:
    """
    Read an NCBI GFF file into a polars dataframe.
    """
    columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]
    return pl.read_csv(
        file,
        has_header=False,
        new_columns=columns,
        separator="\t",
        comment_prefix="#"
    )


def load_seqs(fasta) -> dict[str, str]:
    """
    Load sequences from a fasta file.
    """
    seqs = {}

    with open(fasta) as f:
        current_gene = None
        current_seq = ""
        for line in f:
            if line.startswith(">"):
                if current_gene is not None:
                    seqs[current_gene] = current_seq
                current_gene = line.strip().split()[0].removeprefix(">")
                current_seq = ""
            else:
                current_seq += line.strip()
        
        if current_gene is not None:
            seqs[current_gene] = current_seq
    
    return seqs


def download_genes(model: dict, threads: Optional[int] = None, show_progress: bool = True) -> dict[str, dict]:
    """
    Helper function to download the gene details for all genes in a BiGG model.
    """
    gene_ids = [g["id"] for g in model["genes"]]
    with requests.Session() as s:
        genes = parallel.apply(
            bigg.get_gene_details,
            itertools.repeat(model["id"]),
            gene_ids,
            itertools.repeat(s),
            threads=threads,
            show_progress=show_progress
        )
    return dict(zip(gene_ids, genes))


def download_ncbi_genome(accession: str) -> tuple[dict[str, str], pl.DataFrame]:
    """
    Download protein sequences and genome features from NCBI.

    Returns the sequences and features as dataframes.
    """
    ncbi = NCBI()
    with TemporaryDirectory() as tempdir:
        
        # Download the genome
        download_path = os.path.join(tempdir, "ncbi_data.zip")
        datapkg = ncbi.genome.download(
            accession,
            download_path,
            include_annotation_type=[AnnotationType.PROT_FASTA, AnnotationType.GENOME_GFF],
            validate=True
        )

        with datapkg.zipfile() as archive:

            # Extract the fasta file and load the seqs
            fasta = archive.extract(
                os.path.join("ncbi_dataset", "data", accession, "protein.faa"),
                tempdir
            )
            seqs = load_seqs(fasta)

            # Extract the GFF file and load into dataframe
            gff = archive.extract(
                os.path.join("ncbi_dataset", "data", accession, "genomic.gff"),
                tempdir
            )
            gff_df = read_gff(gff)

    return seqs, gff_df


def format_seq(seq: str, line_length: int = 80):
    """
    Format an amino acid sequence by splitting across multiple lines if longer
    than the specified line length.
    
    Default line length is 80 for writing to a fasta file.
    """
    lines = [seq[i:i+line_length] for i in range(0, len(seq), line_length)]
    return "\n".join(lines)


# Download the model json
model = bigg.get_model(args.model_id)

# Download the gene details
genes = download_genes(model, threads=args.threads, show_progress=True)

# Pull out protein sequences and note any genes that are missing a sequence
missing = []
seqs = {}
for gene_id in genes:
    seq = genes[gene_id].get("protein_sequence")
    if seq is None:
        missing.append(gene_id)
    else:
        seqs[gene_id] = seq
print(f"Retrieved {len(seqs)}/{len(genes)} sequences from BiGG.")

# Try to retreive any missing sequences from NCBI
if len(missing) > 0 and args.refseq_accession is not None:
    ncbi_seqs, ncbi_features = download_ncbi_genome(args.refseq_accession)

    cds = ncbi_features.filter(pl.col.type == "CDS")

    bigg_ncbi_map = (
        pl.DataFrame(
            [genes[x] for x in missing],
            schema_overrides={
                "leftpos": pl.Int64,
                "rightpos": pl.Int64
            }
        )
        .join(
            cds,
            left_on=["genome_name", "leftpos", "rightpos"],
            right_on=["seqid", "start", "end"],
            suffix="_ncbi",
            validate="1:1"
        )
        .with_columns(
            ncbi_id=pl.col.attributes.str.extract(r"ID=cds-(WP_\d+\.\d+)")
        )
        .select("bigg_id", "ncbi_id")
    )
    
    i = 0
    for bigg_id, ncbi_id in bigg_ncbi_map.iter_rows():
        seq = ncbi_seqs.get(ncbi_id)
        if seq is not None:
            seqs[bigg_id] = seq
            i += 1
    print(f"Retrieved {i}/{len(missing)} missing sequences from NCBI.")


# Save protein sequences to output fasta
out = f"{args.model_id}.faa" if args.out is None else args.out
with open(out, "w") as f:
    for gene_id, seq in seqs.items():
        if seq is None:
            continue
        f.write(f">{args.model_id}.{gene_id}\n")
        f.write(f"{format_seq(seq)}\n")
