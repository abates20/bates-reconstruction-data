from typing import Collection, Union, Optional
import warnings

import cobra
from optlang.symbolics import Zero


class pFBAGapfiller:

    def __init__(self, universal: cobra.Model):
        self.universal = universal
    
    def gapfill(
        self,
        model: cobra.Model,
        medium: Optional[dict[str, Union[int, float]]],
        min_fraction: float,
        max_fraction: float,
        constraints: Optional[dict[str, tuple]] = None,
        active_cutoff: float = 1e-6
    ) -> set[str]:
        """
        pFBA gapfiller.

        Modifies universal reaction bag, removes overlapping reacitons from
        universal reaction bag and resets objective if needed, adds model reaction
        to universal bag, sets lower bound for metabolic tasks, sets minimum lower
        bound for previous objective, assemble forward and reverse components of all
        reactions, create objective, based on pFBA, run FBA and identify reactions
        from universal that are now active.
        """

        print("\r[                                         ]", end="", flush=True)

        # Modify universal reaction bag
        with self.universal as universal:  # set the reaction bag as the universal reaction database
            # Set dummy objective to avoid optlang container errors
            universal.objective = Zero

            # Clear any prior medium on universal model
            universal.medium = {}
            
            # Original reaction IDs in model
            orig_rxn_ids = set(r.id for r in model.reactions)

            # Add reactions to universal model (replacing any ones that were already there)
            universal.remove_reactions([r.id for r in model.reactions if universal.reactions.has_id(r.id)])
            universal.add_reactions(r.copy() for r in model.reactions)

            print("\r[---------------                          ]", end="", flush=True)

            # Set medium for model
            if medium is not None:
                universal.medium = medium
            
            # Add provided constraints (reaction bounds)
            if constraints is not None:
                for rxn, bounds in constraints.items():
                    universal.reactions.get_by_id(rxn).bounds = bounds
            
            # Set minimum lower bound for previous objective and add pFBA objective
            universal.objective = model.objective
            pFBAGapfiller._add_pfba(universal, orig_rxn_ids, min_fraction, max_fraction)

            print("\r[-------------------------------          ]", end="", flush=True)

            # Run pFBA and identify reactions from universal that are now active
            solution = universal.optimize()

        active_rxns = solution.fluxes[solution.fluxes.abs() > active_cutoff]
        new_rxn_ids = set(active_rxns.index).difference(orig_rxn_ids)
        print("\r[-----------------------------------------]")

        return new_rxn_ids
    
    def add_reactions(self, model: cobra.Model, new_rxn_ids: Collection[str]) -> cobra.Model:
        """
        Adds new reactions to model by getting reactions and metabolites to be added
        to the model, creates gapfilled model, and identifies extracellular
        metabolites that still need exchanges.
        """

        def get_rxn(rxn_id: str):
            return self.universal.reactions.get_by_id(rxn_id).copy()

        # Add the reactions to the model from the universal model
        model.add_reactions(get_rxn(x) for x in new_rxn_ids)

        return model
    
    @staticmethod
    def _add_pfba(model: cobra.Model, ignore_rxns: Collection[str], min_frac: float, max_frac: float):
        """
        Adds a pFBA objective to the model, but where the objective only factors
        in the flux for reactions not in `ignore_rxns`.
        """
        pFBAGapfiller._set_objective_constraint(model, min_frac, max_frac)

        # Assemble forward and reverse components of all reactions
        coefficientDict = {}
        for rxn in model.reactions:
            rxn: cobra.Reaction
            coeff = int(rxn.id not in ignore_rxns)
            coefficientDict[rxn.forward_variable] = coeff
            coefficientDict[rxn.reverse_variable] = coeff

        # Create objective, based on pFBA
        model.objective = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective.set_linear_coefficients(coefficientDict)

    @staticmethod
    def _set_objective_constraint(model: cobra.Model, min_frac: float, max_frac: float):
        obj_val = max(model.slim_optimize(error_value=None), 0) # ensure objective value isn't < 0
        if obj_val < 1e-6:
            warnings.warn(f"Provided medium for gap-filling does not support biomass production (objective flux = {obj_val:3e})")
        ub, lb = obj_val*max_frac, obj_val*min_frac
        constraint = model.problem.Constraint(model.objective.expression, ub=ub, lb=lb)
        model.add_cons_vars(constraint)

