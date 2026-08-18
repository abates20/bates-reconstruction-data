from pathlib import Path

import cobra

from ..utils import download, DownloadProgress
from ..resources.files import get_resource_path as get_resource_path
from ..resources.managers import (
    universal_manager,
    template_manager,
    diamond_db_manager,
    gene_rxn_map_manager,
    gene_name_manager
)


def get_universal_model(name: str = "universal_prokaryote") -> cobra.Model:
    """
    Get the universal reaction model.

    The first time this function is called, the universal model is loaded from
    the resources directory (which can take 30+ seconds) and then is cached so
    that subsequent calls can simply return the model without loading it again.
    """
    return universal_manager.load(name)


def get_gene_name_map() -> dict[str, str]:
    """
    Get the dictionary mapping KEGG gene IDs to gene names.
    """
    return gene_name_manager.load("kegg")


def get_gene_rxn_map(name: str = "kegg_seed") -> dict[str, list[str]]:
    """
    Get the dictionary mapping KEGG gene IDs to ModelSEED reaction IDs.
    """
    return gene_rxn_map_manager.load(name)


def get_diamond_db_path(name: str = "kegg_prokaryote") -> Path:
    """
    Get the filepath to the KEGG peptide DIAMOND database for blasting.
    """
    return diamond_db_manager.load(name)


def download_diamond_db() -> Path:
    """
    Downloads the DIAMOND database file from the Reconstructor releases page.
    """
    url = "https://github.com/emmamglass/reconstructor/releases/download/v0.0.1/screened_kegg_prokaryotes_pep_db.dmnd"
    progress = DownloadProgress("Downloading the DIAMOND database for blasting...")
    path = get_diamond_db_path("kegg_prokaryote")
    return download(url, path=path, callback=progress)


def remove_diamond_db():
    """
    Deletes the DIAMOND database file. This is meant to be used prior to
    uninstalling Reconstructor; otherwise, the DIAMOND database will be left
    behind when Reconstructor is uninstalled.
    """
    path = get_diamond_db_path()
    path.unlink(missing_ok=True)


def get_template_model(template: str):
    """
    Load one of the bacteria template models.
    """
    return template_manager.load(template)
