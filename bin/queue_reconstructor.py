"""
queue_reconstructor.py
----------------------

Queue a slurm job array to build models with Reconstructor.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("input", type=Path, help="Path to csv listing genomes to be used")
parser.add_argument("--genomes", required=True, type=Path, help="Path to input genomes")
parser.add_argument("--reconstruction", required=True, type=Path, help="Path to reconstruction directory")
parser.add_argument("--models", required=True, type=Path, help="Path to output models directory")
parser.add_argument("--logs", required=True, type=Path, help="Path to logs directory")
parser.add_argument("--gene-rxn-map", required=True, type=Path, help="Path to gene -> reaction map")
parser.add_argument("--diamond-db", required=True, type=Path, help="Path to DIAMOND database")
parser.add_argument("--diamond", help="Path to DIAMOND executable")

media_group = parser.add_argument_group("Media")
media_group.add_argument("--media-db", required=True, type=Path, help="Path to media database tsv")
media_group.add_argument("--gapfill", required=True, help="Media name for gap-filling")
media_group.add_argument("--base-inputs", required=True, help="Media name for base inputs gap-filling")

universe_group = parser.add_argument_group("Universal models")
universe_group.add_argument("--generic-universe", required=True, type=Path, help="Path to generic universal model")
universe_group.add_argument("--gramneg-universe", required=True, type=Path, help="Path to gram-negative universal model")
universe_group.add_argument("--grampos-universe", required=True, type=Path, help="Path to gram-positive universal model")

template_group = parser.add_argument_group("Template models")
template_group.add_argument("--generic-template", required=True, type=Path, help="Path to generic bacteria template model")
template_group.add_argument("--gramneg-template", required=True, type=Path, help="Path to gram-negative template model")
template_group.add_argument("--grampos-template", required=True, type=Path, help="Path to gram-positive template model")

parser.add_argument("--reconstructor-job", required=True, type=Path, help="Path to reconstructor job file")
parser.add_argument("--array-job", required=True, type=Path, help="Path to slurm array job file")
parser.add_argument("--save-jobid", default=None, type=Path, help="Path to file to write slurm job ID")

args = parser.parse_args()


# Additional imports
from typing import Union, Sequence
import subprocess
import warnings
from datetime import datetime
import sys

import polars as pl

from bates_reconstruction_data.taxonomy import Taxonomy, gram_stain


def get_media_cpds(name: str, db: pl.DataFrame) -> list[str]:
    return list(
        db.filter(pl.col.medium == name).get_column("compound") + "_e"
    )

def get_universe_and_template(gramstain: str) -> tuple:
    if gramstain == "negative":
        template = args.gramneg_template
        universe = args.gramneg_universe
    elif gramstain == "positive":
        template = args.grampos_template
        universe = args.grampos_universe
    else:
        template = args.generic_template
        universe = args.generic_universe

    return universe, template


def sbatch(args: Sequence[Union[str, Path]], **sp_kwargs) -> subprocess.CompletedProcess:
    """
    Function for submitting arbitrary slurm jobs.
    """
    cmd = ["sbatch"] + list(args)
    return subprocess.run(cmd, **sp_kwargs)

def fmt(arg) -> str:
    """
    Format a string command line argument and put in quotes.
    """
    return repr(str(arg))


# Load the model info
model_info = pl.read_csv(args.input)

# Add gram stain info for each bacteria
tax_data = {
    "species": [],
    "gram_stain": []
}
for species in model_info.get_column("species").unique():
    tax = Taxonomy.from_name(species)
    gramstain = gram_stain(tax)

    tax_data["species"].append(species)
    if gramstain is None:
        warnings.warn(f"Unable to determine gramstain for {species}. Will use generic template "
        f"{args.generic_template} and universe {args.generic_universe}.")
        tax_data["gram_stain"].append(None)
    else:
        tax_data["gram_stain"].append(gramstain.value)

model_info = (
    model_info.join(
        pl.DataFrame(tax_data),
        on="species",
        how="left"
    )
    .select("organism_id", "refseq_accession", "gram_stain")
)

# Load medium for gap-filling
medium_df = pl.read_csv(args.media_db, separator="\t")
gapfill_cpds = get_media_cpds(args.gapfill, medium_df)
base_inputs = get_media_cpds(args.base_inputs, medium_df)

# Prep folders for storing files
for folder in [args.reconstruction, args.models, args.logs]:
    if not folder.exists():
        folder.mkdir(parents=True)

# Prep command for each bacteria strain
n_jobs = 0
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
command_file = args.reconstruction / f"commands_{timestamp}.txt"

with open(command_file, "w") as f:
    for org_id, accession, gramstain in model_info.iter_rows():
        genome_path = args.genomes / f"{accession}.faa"
        output_path = args.models / f"{accession}.sbml"

        # Check if model has already been built
        if output_path.exists():
            continue

        # Select the template model and universal model based on gram stain
        universe, template = get_universe_and_template(gramstain)

        # Save command to config file
        cmd_args = [
            "bash",
            fmt(args.reconstructor_job),
            "--input_file", fmt(genome_path),
            "--out", fmt(output_path),
            "--media", fmt(",".join(gapfill_cpds)),
            "--base-inputs", fmt(",".join(base_inputs)),
            "--template", fmt(template),
            "--name", fmt(accession.replace(".", "_")),
            "--universal", fmt(universe),
            "--diamond-db", fmt(args.diamond_db),
            "--gene-rxn-map", fmt(args.gene_rxn_map),
            "--blast-output", fmt(args.reconstruction / f"{accession}.blast.out")
        ]

        if args.diamond is not None:
            cmd_args.extend(["--diamond", fmt(args.diamond)])

        f.write(" ".join(cmd_args))
        f.write("\n")
        n_jobs += 1


# Submit job array (if any jobs need to be queued)
if n_jobs > 0:
    log_path = args.logs / f"reconstructor-%A_%a.out"
    slurm_cmd = [
        f"--export=ALL,PYTHON={sys.executable}",
        f"--output={log_path}",
        "--time=1:00:00",
        "--nodes=1",
        "--cpus-per-task=4",
        "--mem=16G",
        f"--array=1-{n_jobs}",
        "--parsable",
        args.array_job,
        command_file
    ]
    res = sbatch(slurm_cmd, stdout=subprocess.PIPE, text=True)
    res.check_returncode()

    # Write job ID to file if specified
    if args.save_jobid is not None:
        with open(args.save_jobid, "w") as f:
            f.write(res.stdout)
