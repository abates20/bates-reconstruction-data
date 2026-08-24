"""
queue_carveme.py
----------------

Queue a slurm job array to build models with CarveMe.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("input", type=Path, help="Path to csv listing genomes to be used")
parser.add_argument("--genomes", required=True, type=Path, help="Path to input genomes")
parser.add_argument("--reconstruction", required=True, type=Path, help="Path to reconstruction directory")
parser.add_argument("--models", required=True, type=Path, help="Path to output models directory")
parser.add_argument("--logs", required=True, type=Path, help="Path to logs directory")
parser.add_argument("--venv", required=True, type=Path, help="Path to script that activates project venv")
parser.add_argument("--gprs", required=True, type=Path, help="Path to CarveMe GPR table")
parser.add_argument("--annotations", required=True, type=Path, help="Path to CarveMe annotation file")
parser.add_argument("--diamond-db", required=True, type=Path, help="Path to DIAMOND database")

media_group = parser.add_argument_group("Media")
media_group.add_argument("--media-db", required=True, type=Path, help="Path to media database tsv")
media_group.add_argument("--gapfill", required=True, help="Media name for gap-filling")

universe_group = parser.add_argument_group("Universal models")
universe_group.add_argument("--generic-universe", required=True, type=Path, help="Path to generic universal model")
universe_group.add_argument("--gramneg-universe", required=True, type=Path, help="Path to gram-negative universal model")
universe_group.add_argument("--grampos-universe", required=True, type=Path, help="Path to gram-positive universal model")

parser.add_argument("--carveme-job", required=True, type=Path, help="Path to carveme job file")
parser.add_argument("--array-job", required=True, type=Path, help="Path to slurm array job file")
parser.add_argument("--save-jobid", default=None, type=Path, help="Path to file to write slurm job ID")

args = parser.parse_args()


# Additional imports
from typing import Union, Sequence
import subprocess
import warnings
from datetime import datetime

import polars as pl

from bates_reconstruction_data.taxonomy import Taxonomy, gram_stain


def get_universe(gramstain: str):
    if gramstain == "negative":
        return args.gramneg_universe
    elif gramstain == "positive":
        return args.grampos_universe
    return args.generic_universe


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
        warnings.warn(f"Unable to determine gramstain for {species}. Will use generic universe "
        f"{args.generic_universe}.")
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

        # Select the universal model based on gram stain
        universe = get_universe(gramstain)

        # Save command to config file
        cmd_args = [
            "bash",
            fmt(args.carveme_job),
            fmt(genome_path),
            "-o", fmt(output_path),
            "--verbose",
            "--universe-file", fmt(universe),
            "--diamond-db", fmt(args.diamond_db),
            "--gpr-file", fmt(args.gprs),
            "--annotation-file", fmt(args.annotations),
            "--blast-output", fmt(args.reconstruction / f"{accession}.blast.out"),
            "--gapfill", fmt(args.gapfill),
            "--mediadb", fmt(args.media_db)
        ]

        f.write(" ".join(cmd_args))
        f.write("\n")
        n_jobs += 1


# Submit job array (if any jobs need to be queued)
if n_jobs > 0:
    log_path = args.logs / f"carveme-%A_%a.out"
    slurm_cmd = [
        f"--export=ENV={args.venv}",
        f"--output={log_path}",
        "--time=1:00:00",
        "--nodes=1",
        "--cpus-per-task=4",
        "--mem=24G",
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
