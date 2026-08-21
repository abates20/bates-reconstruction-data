"""
get_reconstructor_universe.py
-----------------------------

Script to get a universal bacteria model from Reconstructor and correct its
compartments and reaction IDs.
"""

from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("--gramneg", required=True, type=Path, help="Path to output universal gram negative model")
parser.add_argument("--grampos", required=True, type=Path, help="Path to output universal gram positive model")
args = parser.parse_args()


# Additional imports
import cobra

from bates_reconstruction_data.reconstructor import resources
from bates_reconstruction_data.utils import progress


COMPARTMENTS = {
    "c": "cytosol",
    "p": "periplasm",
    "e": "extracellular"
}


# Load universal model
universe = resources.get_universal_model()

# Take '_c' off of reaction IDs
for rxn in progress.track(universe.reactions):
    if rxn.id.endswith("_c"):
        rxn.id = rxn.id.removesuffix("_c")
universe.repair()

# Fix compartments in universal model
universe._compartments = COMPARTMENTS

for met in universe.metabolites:
    if met.compartment == "cytosol":
        met.compartment = "c"
    elif met.compartment == "extracellular":
        met.compartment = "e"
universe.repair()

# Save the gram negative model
with universe as u:
    u.objective = "biomass_GmNeg"
    cobra.io.write_sbml_model(u, args.gramneg)

# Save the gram positive model
with universe as u:
    u.objective = "biomass_GmPos"
    cobra.io.write_sbml_model(u, args.grampos)
