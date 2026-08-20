"""
get_universe.py
-----------------

Script to get a universal bacteria model from CarveMe and correct its
compartments to work with Reconstructor.
"""

from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("url", help="URL to download universal model")
parser.add_argument("-o", "--out", required=True, type=Path, help="Path to output universal model")
args = parser.parse_args()


# Additional imports
from tempfile import TemporaryDirectory

import cobra

from bates_reconstruction_data.utils import download


COMPARTMENTS = {
    "c": "cytosol",
    "p": "periplasm",
    "e": "extracellular"
}


# Download the universal bacteria model from CarveMe
with TemporaryDirectory() as tmp:
    tmpf = Path(tmp) / "universe_bacteria.xml.gz"
    download.download(args.url, tmpf)
    universe = cobra.io.read_sbml_model(tmpf)

# Fix the compartments in the universal model
universe._compartments = COMPARTMENTS

for met in universe.metabolites:
    met: cobra.Metabolite
    if met.compartment is not None and met.compartment.startswith("C_"):
        met.compartment = met.compartment.removeprefix("C_")

# Save the corrected model
cobra.io.write_sbml_model(universe, args.out)
