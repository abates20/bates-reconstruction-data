from typing import Union, Literal, Sequence, Optional
import warnings
import os
import enum
from pathlib import Path

import cobra
from cobra.io.sbml import CobraSBMLError

from . import resources, _funcs, utils
from .gapfill import pFBAGapfiller
from .medium import get_medium


class FileType(enum.Enum):
    FASTA = 1
    BLAST_OUTPUT = 2
    SBML = 3


def reconstruct(
    input_file: Union[str, os.PathLike],
    *,
    file_type: Literal[1, 2, 3] = 1,
    media: Union[Literal["rich", "minimal", "complete"], Sequence[str]] = "rich",
    base_inputs: Optional[Union[Literal["rich", "minimal", "complete"], Sequence[str]]] = None,
    org: Optional[str] = None,
    min_frac: float = 0.01,
    max_frac: float = 1.0,
    template: Union[cobra.Model, Literal["gram_negative", "gram_positive"]] = "gram_negative",
    model_id: Optional[str] = None,
    universal: Union[cobra.Model, str] = "universal_prokaryote",
    diamond_db: Union[str, os.PathLike] = "kegg_prokaryote",
    gene_rxn_map: Union[dict[str, list[str]], str] = "kegg_seed",
    cpu: Optional[int] = None,
    gapfill: bool = True,
    open_exchanges: bool = True,
    blast_output: Optional[Union[str, bytes, os.PathLike]] = None,
    diamond_path: Optional[Union[str, bytes, os.PathLike]] = None
) -> cobra.Model:
    """
    Reconstruct a genome-scale metabolic model from an input file.

    For a FASTA input, uses DIAMOND to find similar protein sequences in KEGG.
    The blast output is then used to create a draft metabolic model, which is
    gap-filled with a pFBA-based approach.

    Parameters
    ----------
    input_file : str | PathLike
        The file path to the input file. The input file should be 1) a FASTA
        with amino acid sequences, 2) a tsv output from DIAMOND blast, or 3) an
        SBML model file.

    file_type : Literal[1, 2, 3]
        1 for a FASTA file, 2 for a DIAMOND blast output file, or 3 for an SBML
        model file. The default is 1.

    media : str | Sequence[str]
        Either the name of a registered media composition ("rich", "minimal", or
        "complete") or a sequence of extracellular metabolite IDs. The default
        is "rich".

    base_inputs : str | Sequence[str], optional
        Either the name of a registered media composition ("rich", "minimal", or
        "complete") or a sequence of extracellular metabolite IDs. This media is
        used for a second round of gap-filling where the model is forced to
        uptake each of the metabolites in the media. If `None`, the second round
        of gap-filling is skipped. Default is `None`.

    org : str, optional
        A KEGG organism ID. If supplied, all genes from this organism will added
        to the DIAMOND blast results. Default is None.

    min_frac : float
        The minimum objective fraction to use for pFBA gap-filling. Default is
        0.01.

    max_frac : float
        The maximum objective fraction to use for pFBA gap-filling. Default is
        0.5.

    template : cobra.Model | str
        The template model to use when creating a new model. Can either be a
        `cobra.Model` object or the gram stain of the bacteria for which a model
        is being reconstructed. If gram stain is given, the corresponding model
        template will be selected. Default is "gram_negative".

    model_id : str, optional
        The ID to use for the model that is created.

    universal : str | cobra.Model
        The universal model that reactions will be pulled from and that will be
        used for gap-filling. Should either be a cobra.Model object or the name
        of a registered universal model.

    diamond_db : str | os.PathLike
        The diamond database for blasting. Should either be a path to a database
        file or the name of of a registered database.

    gene_rxn_map : str | dict
        A dictionary mapping gene IDs to reaction IDs. Should either be a
        dictionary object or the name of a registered map.

    cpu : int, optional
        How many processors to use for DIAMOND. Default is None, which means
        that DIAMOND will automatically detect how many to use.

    gapfill : bool
        Whether or not to perform pFBA gap-filling on the model. Default is
        True.

    open_exchanges : bool
        Whether or not to set the exchanges of the returned model to open or
        not. Default is True.

    blast_output : str | os.PathLike, optional
        A path to where the DIAMOND blast output should be saved if the input
        file is an annotated genome. Defaults to the path of the annotated
        genome with the stuffix '.KEGGProt.out'.

    diamond_path : str | os.PathLike, optional
        A path to the DIAMOND executable to use.

    Returns
    -------
    cobra.Model
        The reconstructed metabolic model as a COBRApy model.
    """

    # Convert input_file to path and ensure it exists
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(
            f"The provided input file path does not exist: {input_file}"
        )

    # Convert file type to enum value
    ftype = FileType(file_type)

    # Ensure min/max fractions are properly set
    if min_frac <= 0.0 or min_frac > 1.0:
        warnings.warn(
            f"Improper minimum fraction selected ({repr(min_frac)}). Defaulting to 0.01"
        )
        min_frac = 0.01
    if max_frac <= 0.0 or max_frac > 1.0:
        warnings.warn(
            f"Improper maximum fraction selected ({repr(max_frac)}). Defaulting to 1.0"
        )
        max_frac = 1.0
    if max_frac < min_frac:
        warnings.warn(
            f"Input maximum fraction ({repr(max_frac)}) less than minimum fraction ({repr(min_frac)}). Minimum set to half maximum"
        )
        min_frac = max_frac * 0.5

    # Determine the media for gap-filling
    if isinstance(media, str):
        media = get_medium(media)

    # Load universal reaction model
    if not isinstance(universal, cobra.Model):
        print("Loading universal model...")
        universal = resources.get_universal_model(universal)

    if ftype in [FileType.FASTA, FileType.BLAST_OUTPUT]:

        # Run BLAST if needed
        if ftype == FileType.FASTA:
            if not os.path.exists(diamond_db):
                diamond_db = resources.get_diamond_db_path(diamond_db)
            print(
                f"Blasting input sequences against {diamond_db}, may take some time..."
            )
            if blast_output is None:
                blast_output = input_path.with_suffix(".KEGGprot.out")
            _funcs.run_blast(
                input_path,
                blast_output,
                diamond_db,
                cpu,
                diamond_path
            )

        else:
            blast_output = input_path

        # Draft model from BLAST output
        print("Drafting model from BLAST results...")
        gene_hits = _funcs.read_blast(blast_output)

        # Load the gene -> rxn map if needed
        if isinstance(gene_rxn_map, str):
            gene_rxn_map = resources.get_gene_rxn_map(gene_rxn_map)

        # Add any extra genes from a provided organism
        if org is not None:
            gene_hits = _funcs.add_genes_from_organism(gene_hits, org, gene_rxn_map)

        # Create the rxn -> gene map
        rxns = _funcs.genes_to_rxns(gene_hits, gene_rxn_map)

        # Create gene objects
        gene_names = resources.get_gene_name_map()
        genes = _funcs.create_gene_objects(gene_hits, gene_names)

        # Draft the model
        draft_model = _funcs.create_model(template, rxns, genes, universal, model_id)

    else:
        try:
            draft_model = cobra.io.read_sbml_model(input_path)
        except CobraSBMLError:
            draft_model = cobra.io.load_json_model(input_path)

    # Reactions/metabolites in draft model
    draft_reactions = set([x.id for x in draft_model.reactions])
    draft_metabolites = set([x.id for x in draft_model.metabolites])

    if gapfill:
        print("Gapfilling model...")
        gapfiller = pFBAGapfiller(universal)

        # Deteremine the medium for gap-filling
        medium_dict = {}
        for cpd in media:
            exch_id = f"EX_{cpd}"
            if universal.reactions.has_id(exch_id) or draft_model.reactions.has_id(exch_id):
                medium_dict[exch_id] = 1000.0

        # Gap-filling
        new_reactions = gapfiller.gapfill(draft_model, medium_dict, min_frac, max_frac)
        filled_model = gapfiller.add_reactions(draft_model, new_reactions)

        # 2nd round of gap-filling if base_inputs provided
        if base_inputs is not None:
            if isinstance(base_inputs, str):
                base_inputs = get_medium(base_inputs)

            constraints = {f"EX_{cpd}": (-1000.0, -0.01) for cpd in base_inputs}
            media_reactions = gapfiller.gapfill(
                filled_model,
                medium_dict,
                min_frac,
                max_frac,
                constraints
            )
            filled_model = gapfiller.add_reactions(filled_model, media_reactions)

        # Fill in annotations
        final_model = _funcs.add_annotation(filled_model)

    else:
        final_model = _funcs.add_annotation(draft_model)

    # Make sure all extracellular metabolites have exchanges
    utils.add_missing_exchanges(final_model)

    if not open_exchanges:
        for exch in final_model.exchanges:
            exch.bounds = (0.0, 0.0)
    elif ftype != FileType.SBML:
        media_exchanges = set(f"EX_{cpd}" for cpd in media)
        for exch in final_model.exchanges:
            if exch.id in media_exchanges:
                exch.bounds = (-1000.0, 1000.0)
            else:
                exch.bounds = (0.0, 1000.0)

    # Check the final model
    _funcs.check_model(draft_reactions, draft_metabolites, final_model)

    return final_model
