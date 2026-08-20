"""
make_bigg_template.py
----------------------

Make template model for the BiGG namespace (based on CarveMe).
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("biomass_db", type=Path, help="Path to biomass_db.tsv specifying biomass compositions")
parser.add_argument("-o", "--out", type=Path, required=True, help="Path to save template model")
parser.add_argument("--biomass", required=True, help="Which biomass formulation to use from biomass db")
parser.add_argument("--template-id", default=None, help="Model ID for template model")
parser.add_argument("--threads", type=int, default=None, help="Number of threads for parallel execution")
args = parser.parse_args()


# Additional imports
import re

import polars as pl
import cobra

from bates_reconstruction_data import bigg
from bates_reconstruction_data.utils import parallel


COMPARTMENTS = {
    "c": "cytosol",
    "p": "periplasm",
    "e": "extracellular"
}


def create_metabolite(bigg_id, compartment) -> cobra.Metabolite:
    """
    Create a metabolite object by retrieving its info from BiGG.
    """

    # Get info about metabolite from BiGG
    met_info: dict = bigg.get_universal_metabolite(bigg_id)

    formulae = met_info.get("formulae", [])
    formula = formulae[0] if len(formulae) > 0 else None

    charges = met_info.get("charges", [])
    charge = charges[0] if len(charges) > 0 else None

    # Create metabolite
    met = cobra.Metabolite(
        f"{bigg_id}_{compartment}",
        formula,
        met_info.get("name"),
        charge,
        compartment
    )
    add_identifiers(met, met_info.get("database_links", {}))
    
    return met


def add_identifiers(met: cobra.Metabolite, db_dict: dict[str, list[dict[str, str]]]) -> None:
    """
    Add links for other databases to a metabolite's annotation dict.
    """
    for identifiers in db_dict.values():
        key = None
        for identifier in identifiers:
            id_ = identifier.get("id")

            if key is None:
                keys = re.findall(
                    f"^http://identifiers.org/(.+)/{id_}$",
                    identifier.get("link", "")
                )
                if len(keys) > 0:
                    key = keys[0] # There should only be one group from regex
                    met.annotation[key] = []
                else:
                    continue

            met.annotation[key].append(id_)


# Create template model
template = cobra.Model(args.template_id)
template.compartments = COMPARTMENTS

# Load biomass data
biomass_db = (
    pl.read_csv(args.biomass_db, separator="\t")
    .select("bigg_id", "comp", args.biomass)
    .drop_nulls()
)

# Create metabolites for biomass reaction
met_list = parallel.apply(
    create_metabolite,
    biomass_db.get_column("bigg_id"),
    biomass_db.get_column("comp"),
    threads=args.threads
)
metabolites = {m.id: m for m in met_list} # turn into dict for easy lookup

# Prep stoichiometry for biomass reaction
coeffs = {}
for bigg_id, compartment, coeff in biomass_db.iter_rows():
    met = metabolites[f"{bigg_id}_{compartment}"]
    coeffs[met] = coeff

# Create biomass reaction and add to model
biomass = cobra.Reaction("biomass", name="Biomass reaction", lower_bound=0, upper_bound=1000)
biomass.add_metabolites(coeffs)
template.add_reactions([biomass])
template.objective = biomass

# Save template
cobra.io.write_sbml_model(template, args.out)
