"""
get_reconstructor_templates.py
------------------------------

Script to prepare gram positive and gram negative template models based on the
ModelSEED universal model.
"""


from argparse import ArgumentParser
from pathlib import Path


parser = ArgumentParser()
parser.add_argument("universe", help="Path to ModelSEED universal model")
parser.add_argument("--gramneg", required=True, type=Path, help="Path to output gram negative template")
parser.add_argument("--grampos", required=True, type=Path, help="Path to output gram positive template")
args = parser.parse_args()


# Additional imports
import cobra


COMPARTMENTS = {
    "c": "cytosol",
    "p": "periplasm",
    "e": "extracellular"
}


universe = cobra.io.read_sbml_model(args.universe)

# Create a template model for gram negative bacteria
gram_neg = cobra.Model("gram_neg", "Gram negative template")

# Specify compartments
gram_neg.compartments = COMPARTMENTS

# Add biomass reactions
gram_neg_rxns = [
    "biomass_GmNeg",
    "EX_biomass",
    "cofactor_rxn",
    "lipid_rxn",
    "protein_rxn",
    "dna_rxn",
    "rna_rxn"
]
for rxn_id in gram_neg_rxns:
    rxn: cobra.Reaction = universe.reactions.get_by_id(rxn_id).copy()
    rxn.annotation["sbo"] = "SBO:0000629"
    gram_neg.add_reactions([rxn])

# Set biomass objective
gram_neg.objective = "biomass_GmNeg"

# Save model
cobra.io.write_sbml_model(gram_neg, args.gramneg)

# Create a template model for gram positive bacteria
gram_pos = cobra.Model("gram_pos", "Gram positive template")

# Specify compartments
gram_pos.compartments = COMPARTMENTS

# Add biomass reactions
gram_pos_rxns = [
    "biomass_GmPos",
    "EX_biomass",
    "cofactor_rxn",
    "lipid_rxn",
    "protein_rxn",
    "dna_rxn",
    "rna_rxn",
    "GmPos_cellwall", # doesn't show up as mandatory precursor because of GmNeg_cellwall reaction
    "teichoicacid_rxn",
    "peptidoglycan_rxn"
    
]
for rxn_id in gram_pos_rxns:
    rxn: cobra.Reaction = universe.reactions.get_by_id(rxn_id).copy()
    rxn.annotation["sbo"] = "SBO:0000629"
    gram_pos.add_reactions([rxn])

# Set biomass objective
gram_pos.objective = "biomass_GmPos"

# Save model
cobra.io.write_sbml_model(gram_pos, args.grampos)
