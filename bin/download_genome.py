"""
download_genome.py
------------------

This script downloads the annotated genomes from NCBI that will be used to
build models throughout this project.
"""

from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("accession", help="RefSeq accession to download")
parser.add_argument("-o", "--out", required=True, type=Path, help="Path to output fasta")
args = parser.parse_args()


# More imports
from tempfile import TemporaryDirectory
import shutil
import os

from bates_reconstruction_data.ncbi_genomes import NCBI, AnnotationType


# Download the genome from NCBI and put it in the appropriate directory
ncbi = NCBI()
with TemporaryDirectory() as tempdir:
    
    # Download the genome
    download_path = os.path.join(tempdir, "ncbi_data.zip")
    datapkg = ncbi.genome.download(
        args.accession,
        download_path,
        include_annotation_type=AnnotationType.PROT_FASTA,
        validate=True
    )

    # Put the FASTA file in the appropriate spot
    with datapkg.zipfile() as archive:
        file = os.path.join("ncbi_dataset", "data", args.accession, "protein.faa")
        tmp_path = archive.extract(file, tempdir)
        shutil.move(tmp_path, args.out)
