"""
build_model_set.py
------------------

Build a set of models for an organism using the provided BLAST results and a
table specifying the annotation set to use for each model.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("input", type=Path, help="Path to file specifying annotation sets")
parser.add_argument("--out", type=Path, required=True, help="Path to output directory")
parser.add_argument("--blast", type=Path, required=True, help="Path to blast results file")
parser.add_argument("--universe", type=Path, required=True, help="Path to universal model")
parser.add_argument("--gene-rxn-map", type=Path, required=True, help="Path to base gene -> reaction map")
parser.add_argument("--template", type=Path, required=True, help="Path to template model")
parser.add_argument("--mediadb", type=Path, required=True, help="Path to tsv with media for gap-filling")
parser.add_argument("--processes", type=int, default=None, help="Number of parallel processes")
parser.add_argument("--solver", default=None, help="What solver to use")
parser.add_argument("--log", default="INFO", help="Logging level")
parser.add_argument("--rebuild", action="store_true", help="Rebuild existing models")
args = parser.parse_args()


# Additional imports
import gzip
from itertools import repeat
import json
import logging
from tempfile import TemporaryDirectory
from typing import Sequence, Optional

import polars as pl
import cobra

from bates_reconstruction_data.reconstructor.build import reconstruct
from bates_reconstruction_data.utils import parallel


# Media for gap-filling
MINIMAL_MEDIA = "minimal[bigg]"
BASE_INPUTS = "complete[bigg]"

# Configure logging
logging.basicConfig(
    format="[%(asctime)s] : %(levelname)s : %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=args.log
)

# Set cobra solver
COBRA_CONFIG = cobra.Configuration()
if args.solver is not None:
    COBRA_CONFIG.solver = args.solver
logging.info(f"Cobra solver set to {COBRA_CONFIG.solver}")

# Log the provided inputs
logging.info(f"Input annotation sets: {args.input}")
logging.info(f"Ouptut directory: {args.out}")
logging.info(f"Input BLAST results: {args.blast}")
logging.info(f"Universe file: {args.universe}")
logging.info(f"Base gene -> rxn map: {args.gene_rxn_map}")
logging.info(f"Template model: {args.template}")


def get_media_cpds(name: str, db: pl.DataFrame) -> list[str]:
    """
    Get the set of compounds in a specific medium.
    """
    return list(
        db.filter(pl.col.medium == name).get_column("compound") + "_e"
    )


def open_map(path: Path, mode: str = "rt"):
    """
    Opens a json file (that may be gzipped).
    """
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def build_model(
    blast: Path,
    out: Path,
    universe: Path,
    gene_rxn_map: Path,
    template: Path,
    medium: Sequence[str],
    base_inputs: Sequence[str]
) -> None:
    """
    Build a model from the specified blast file.
    """
    logging.info(f"Building model from input {blast}")
    logging.info(f"Using gene -> rxn map {gene_rxn_map}")

    with open_map(gene_rxn_map) as f:
        map = json.load(f)

    model = reconstruct(
        blast,
        file_type=2,
        media=medium,
        base_inputs=base_inputs,
        template=cobra.io.read_sbml_model(template),
        universal=cobra.io.read_sbml_model(universe),
        gene_rxn_map=map
    )

    logging.info(f"Finished building model from {gene_rxn_map}. Saving to {out}")
    cobra.io.save_json_model(model, out)


# Get compound lists for gap-filling
media_db = pl.read_csv(args.mediadb, separator="\t")
medium = get_media_cpds(MINIMAL_MEDIA, media_db)
base_inputs = get_media_cpds(BASE_INPUTS, media_db)

# Load input annotation sets
annotation_sets = (
    pl.read_csv(
        args.input,
        schema_overrides={"set_id": pl.String}
    )
    .with_columns(
        pl.col.new_function.str.split(";")
    )
)

# Load base gene -> rxn map file
with open_map(args.gene_rxn_map) as f:
    base_gene_rxn_map = json.load(f)

# Get annotation sets that models haven't been built for yet and save
# corresponding gene -> rxn maps to temp location
with TemporaryDirectory() as tmpdir:
    logging.debug(f"tmpdir = {tmpdir}")

    tmp_maps = []
    outputs = []
    for set_id, gene, new_function in annotation_sets.iter_rows():
        out = args.out / f"{set_id}.json"
        map_path = Path(tmpdir) / f"{set_id}.json"

        if args.rebuild or not out.exists():
            tmp_maps.append(map_path)
            outputs.append(out)

            # Create new gene -> rxn map
            if gene is None:
                new_map = base_gene_rxn_map
            else:
                new_map = base_gene_rxn_map.copy()
                new_map[gene] = new_function

            # Save temporary gene -> rxn map file
            with open(map_path, "w") as f:
                json.dump(new_map, f)

    logging.info(f"{len(tmp_maps)} out of {annotation_sets.height} models need to be built")

    # Build models
    parallel.apply(
        build_model,
        repeat(args.blast),
        outputs,
        repeat(args.universe),
        tmp_maps,
        repeat(args.template),
        repeat(medium),
        repeat(base_inputs),
        max_workers=args.processes,
        executor_type="process",
        show_progress=False
    )
