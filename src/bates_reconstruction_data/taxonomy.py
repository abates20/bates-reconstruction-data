"""
taxonomy
--------

A python module for using information from the NCBI Taxonomy database.
"""


from typing import Callable, TypeVar, ParamSpec, Optional, Union, IO, Sequence, Any, Iterable
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import io
from functools import wraps
import warnings
from enum import Enum
import os
import time

import numpy as np
from numpy.typing import NDArray
import polars as pl

from bates_reconstruction_data.utils import download, atomicfile


_TAXONOMY_ORDER = [
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]
_GRAM_POSITIVE_PHYLA = [
    "Actinobacteria",
    "Chloroflexi",
    "Firmicutes",
    "Tenericutes",
    "Actinomycetota",
    "Chloroflexota",
    "Bacillota",
    "Mycoplasmatota",
]

_GRAM_NEGATIVE_CLASS_EXCEPTIONS = ["Negativicutes"]
_RECORD_SEP = "\x1E"
_DMP_COLUMNS = [
    "tax_id",
    "tax_name",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
    "domain"
]

_T = TypeVar("_T")
_P = ParamSpec("_P")


def _cache_single(func: Callable[_P, _T]) -> Callable[_P, _T]:
    """
    Caches the result for a function that always returns the same value
    regardless of the input parameters.
    """
    called = False
    result = None

    @wraps(func)
    def wrapper(*args: _P.args, **kwargs: _P.kwargs) -> _T:
        nonlocal called, result

        if not called:
            result = func(*args, **kwargs)
            called = True
        return result  # ty:ignore[invalid-return-type]

    return wrapper

@_cache_single
def get_taxonomy_data(local: Optional[Union[str, bytes, Path]] = None) -> pl.DataFrame:
    """
    Load the table of taxonomic lineages from NCBI.

    If this table has not already been loaded and cached then it gets downloaded
    from NCBI if needed or loaded from the provided local path.
    """
    if local is not None:
        source = local
        if not os.path.exists(local):
            with atomicfile.AtomicFileWrapper(local, "b") as buffer:
                _download_taxonomy_data(buffer)
    else:
        source = io.BytesIO()
        _download_taxonomy_data(source)
        source.seek(0)

    return pl.read_csv(
        source,
        has_header=False,
        new_columns=_DMP_COLUMNS,
        separator=_RECORD_SEP,
        quote_char=None
    )


def distance_matrix(taxa: Sequence["Taxonomy"]) -> NDArray[np.float64]:
    """
    Create a taxonomic distance matrix for a set of taxa.
    """
    n = len(taxa)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if j < i:
                dist_matrix[i,j] = dist_matrix[j,i]

            elif i == j:
                dist_matrix[i,j] = 0

            else:
                d = taxonomic_distance(taxa[i], taxa[j])
                dist_matrix[i,j] = d

    return dist_matrix


def taxonomic_distance(t1: "Taxonomy", t2: "Taxonomy") -> int:
    """
    Determine the taxonomic distance between two taxa.
    """
    i = 0
    for level in _TAXONOMY_ORDER:
        if t1[level] != t2[level]:
            break
        i += 1
    return len(_TAXONOMY_ORDER) - i


def lowest_equal_taxonomic_level(t1: "Taxonomy", t2: "Taxonomy") -> Optional[str]:
    """
    Find the most granular taxonomic level that is shared between two taxa.
    """
    level = None
    for next_level in _TAXONOMY_ORDER:
        if t1[next_level] != t2[next_level]:
            break
        level = next_level
    return level


def taxonomic_diversity(species: Iterable[str]) -> float:
    """
    Calculate the taxonomic diversity for a set of species.

    Taxonomic diversity is defined here as the average distance between any two
    species.
    """
    species = list(set(species)) # use set to get just unique species
    n = len(species)
    if n <= 1: # Prevent divide by 0 error
        return 0

    # Get a taxonomy object for each species
    taxa = get_taxa(species)

    # Create distance matrix and calculate diversity
    dist_matrix = distance_matrix(taxa)
    summation = 0
    for j in range(n):
        for i in range(j):
            summation += dist_matrix[i,j]
    return summation / (n * (n-1) / 2)


def unweighted_coverage(species: Iterable[str], ref_species: Iterable[str]) -> float:
    """
    Calculate the unweighted coverage for a set of species compared to a
    reference set of species.

    This is simply the fraction of species from the reference set that are
    included in the input set.
    """
    unique_input = set(species)
    unique_ref = set(ref_species)
    num_covered = len(unique_input & unique_ref)
    return num_covered / len(unique_ref)


def weighted_coverage(species: Iterable[str], ref_species: Iterable[str], max_distance: float = 8) -> float:
    """
    Calculate the weighted coverage for a set of species compared to a reference
    set of species.

    Max distance defaults to 8 because that is the maximum possible distance
    between two taxa (i.e., taxa from different domains).
    """
    # Get taxonomy object for all species (fast)
    unique_input = set(species)
    unique_ref = set(ref_species)
    all_species = list(unique_input | unique_ref)
    taxa = get_taxa(all_species)

    # Build distance matrix (slow)
    dist_matrix = distance_matrix(taxa)

    # Subset distance matrix so rows are input species and columns are ref species (fast)
    row_idx = [i for i, sp in enumerate(all_species) if sp in unique_input]
    col_idx = [i for i, sp in enumerate(all_species) if sp in unique_ref]
    dist_matrix = dist_matrix[np.ix_(row_idx, col_idx)]

    # Find minimum distance for each ref species to input species (fast)
    min_dists: NDArray[np.float64] = np.min(dist_matrix, axis=0)

    # Average of minimum distances (fast)
    avg_dist = np.mean(min_dists)

    return float(1 - (avg_dist / max_distance))


def get_taxa(names: Iterable[str], local: Optional[Union[str, bytes, Path]] = None) -> list["Taxonomy"]:
    """
    Create a taxonomy object for each provided taxonomic name.
    """
    tax_data: pl.DataFrame = (
        get_taxonomy_data(local)
        .join(
            pl.DataFrame({"tax_name": names}),
            on=pl.col("tax_name").str.to_lowercase(),
            how="right",
            maintain_order="right" # Keep resulting taxonomy objects in same order as provided names
        )
    )

    # Check for null rows
    nulls = tax_data.filter(pl.col("tax_id").is_null())
    if not nulls.is_empty():
        raise NoTaxonomyDataError(names=list(nulls.get_column("tax_name_right")))

    # Remove tax_name_right column (doesn't coalesce because joining with expression)
    tax_data = tax_data.drop("tax_name_right")

    taxa = [Taxonomy.from_row_dict(row) for row in tax_data.to_dicts()]
    return taxa


def _clean_dmp(dmp: IO[bytes], out: IO[bytes]):
    replacements = [
        ("\t|\t".encode(), _RECORD_SEP.encode()),
        ("\t|\n".encode(), "\n".encode())
    ]
    for line in dmp:
        for old, new in replacements:
            line = line.replace(old, new)
        out.write(line)


def _download_taxonomy_data(buffer: IO[bytes]) -> None:
    url = "https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip"
    with TemporaryDirectory() as tempdir:
        tempdir = Path(tempdir)
        path = tempdir / "taxdump.zip"
        download.download(url, path, download.DownloadProgress("Downloading taxonomy data..."))


        with ZipFile(path) as zf:
            _clean_dmp(zf.open("rankedlineage.dmp"), buffer)


class TaxonomyError(Exception):

    @classmethod
    def enforce_error_type(cls, func):

        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                result = func(*args, **kwargs)
            except Exception as error:
                if not isinstance(error, cls):
                    raise cls() from error
                else:
                    raise
            else:
                return result

        return wrapper


class TaxonomyFileNotFoundError(FileNotFoundError, TaxonomyError):

    def __init__(self):
        msg = "No taxonomy data file was found."
        super().__init__(msg)


class NoTaxonomyDataError(TaxonomyError):

    def __init__(self, **kwargs):
        constraints = [f"{k}={repr(v)}" for k, v in kwargs.items()]
        msg = f"Taxonomy information was not found for the provided input: {', '.join(constraints)}"
        super().__init__(msg)


class DmpParseError(TaxonomyError):

    def __init__(self):
        msg = "An error occured parsing the taxonomy .dmp file from NCBI"
        super().__init__(msg)


class Taxonomy:
    """
    Represents a taxonomic classification from the NCBI
    taxonomy database.

    Attributes
    ----------
    name : str, optional
        The taxonomy name.
    ncbi_id : int, optional
        The NCBI ID of the taxonomy.
    species_ : str, optional
        The species of the taxonomy.
    genus_ : str, optional
        The genus of the taxonomy.
    family_ : str, optional
        The family of the taxonomy.
    order_ : str, optional
        The order of the taxonomy.
    class_ : str, optional
        The class of the taxonomy.
    phylum_ : str, optional
        The phylum of the taxonomy.
    kingdom_ : str, optional
        The kingdom of the taxonomy.
    domain_ : str, optional
        The domain of the taxonomy.
    """

    def __init__(
        self,
        *,
        name: Optional[str] = None,
        ncbi_id: Optional[int] = None,
        species_: Optional[str] = None,
        genus_: Optional[str] = None,
        family_: Optional[str] = None,
        order_: Optional[str] = None,
        class_: Optional[str] = None,
        phylum_: Optional[str] = None,
        kingdom_: Optional[str] = None,
        domain_: Optional[str] = None,
    ):
        self._name: Optional[str] = name
        self.ncbi_id: Optional[int] = ncbi_id
        self.species_: Optional[str] = species_
        self.genus_: Optional[str] = genus_
        self.family_: Optional[str] = family_
        self.order_: Optional[str] = order_
        self.class_: Optional[str] = class_
        self.phylum_: Optional[str] = phylum_
        self.kingdom_: Optional[str] = kingdom_
        self.domain_: Optional[str] = domain_

    def __repr__(self):
        args = {k: v for k, v in vars(self).items() if v is not None}
        name = args.pop("_name", None)

        kwargs = ", ".join([f"{k}={repr(v)}" for k, v in args.items()])

        return (
            f"Taxonomy({name}, {kwargs})" if name is not None else f"Taxonomy({kwargs})"
        )

    def __str__(self):
        elbow = "\u2514"
        hline = "\u2500"
        lines = [f"TaxonomyID: {self.ncbi_id}"]

        # Find most granular level
        for level in reversed(_TAXONOMY_ORDER):
            if level is not None:
                final_level = level
                break

        # Build lines
        for i, level in enumerate(_TAXONOMY_ORDER):
            name = self.__getattribute__(f"{level}_")
            indent = "   " * i
            line = f"{indent}{elbow}{hline} {level.title()}: {name}"
            lines.append(line)
            if level == final_level:
                break

        return "\n".join(lines)

    def __getitem__(self, key: str):
        if key in _TAXONOMY_ORDER:
            key = f"{key}_"
        return self.__getattribute__(key)

    @property
    def name(self):
        """
        The name for the taxonomic classification if provided,
        otherwise the species name.
        """
        return self._name if self._name is not None else self.species_name

    @property
    def species_name(self):
        """
        The species for the taxonomic classification if provided,
        otherwise the genus + 'sp.'.
        """
        if self.species_ is not None:
            return self.species_
        genus = self.genus_ if self.genus_ is not None else "[Unknown genus]"
        return f"{genus} sp."

    @classmethod
    def from_name(cls, name: str) -> "Taxonomy":
        """
        Creates a taxonomy object by searching the NCBI taxonomy database for
        the provided name.
        """
        data: pl.DataFrame = (
            get_taxonomy_data()
            .filter(pl.col("tax_name").str.to_lowercase().eq(name.lower()))
        )

        if data.is_empty():
            raise NoTaxonomyDataError(name=name)

        if data.height > 1:
            raise TaxonomyError("Input query matched multiple rows in taxonomy data")

        taxonomy = data.to_dicts()[0]
        return cls.from_row_dict(taxonomy)

    @classmethod
    def from_id(cls, id_: int) -> "Taxonomy":
        data: pl.DataFrame = (
            get_taxonomy_data()
            .filter(pl.col("tax_id").eq(id_))
        )

        if data.is_empty():
            raise NoTaxonomyDataError(tax_id=id_)

        if data.height > 1:
            raise TaxonomyError("Input query matched multiple rows in taxonomy data")

        taxonomy = data.to_dicts()[0]
        return cls.from_row_dict(taxonomy)

    @classmethod
    def from_row_dict(cls, data: dict[str, Any]) -> "Taxonomy":
        """
        Create a taxonomy object from a row in the NCBI taxonomy database
        (provided as a dictionary).
        """
        # Take the tax_name and put that in as the first null level
        prev = None
        for level in reversed(_TAXONOMY_ORDER):
            if data.get(level) is not None:
                if prev is not None:
                    data[prev] = data.get("tax_name")
                break
            prev = level

        # Change names to what the Taxonomy object expects
        data["name"] = data.pop("tax_name", None)
        data["ncbi_id"] = data.pop("tax_id", None)
        for level in _TAXONOMY_ORDER:
            data[f"{level}_"] = data.pop(level, None)

        return cls(**data)


class Gram(Enum):
    POSITIVE = "positive"
    NEGATIVE = "negative"


def gram_stain(taxonomy: Union[str, Taxonomy]) -> Optional[Gram]:
    """
    Determine whether a bacteria is gram negative or positive.

    Gets taxonomic information from the input information using
    `bactax.get_taxonomy` and then determines the gram stain
    based on the assigned phylum/class. The method was taken
    from the [AMR R package](https://msberends.github.io/AMR/)
    where any bacteria in one of the following phyla are
    considered gram positive:

    - Actinomycetota (previously called Actinobacteria)
    - Chloroflexota (previously called Chloroflexi)
    - Bacillota (previously called Firmicutes) ... except the Negativicutes class!
    - Mycoplasmatota (previously called Tenericutes)

    Since the gram stain is based on phylum/class information,
    at least class (or more granular) must be provided as input to
    avoid `None` being returned.

    Returns
    -------
    Gram.NEGATIVE | Gram.POSITIVE | None
        Returns either `Gram.NEGATIVE` or `Gram.POSITIVE` if the
        phylum/class can be determined. If the phylum/class are not
        found, then it returns `None` (because the gram stain cannot
        be determined).
    """

    if isinstance(taxonomy, str):
        taxonomy = Taxonomy.from_name(taxonomy)

    # Ensure the taxonomy is for a bacteria
    if taxonomy.domain_ != "Bacteria":
        warnings.warn(f"Gram stain does not apply to domain/realm {taxonomy.domain_}.")

    if taxonomy.phylum_ is None or (
        taxonomy.phylum_ == "Bacillota" and taxonomy.class_ is None
    ):
        warnings.warn(f"Unable to determine gram stain for Bacillota phylum with no class information.")
        return None

    is_gram_positive = (
        taxonomy.phylum_ in _GRAM_POSITIVE_PHYLA
        and taxonomy.class_ not in _GRAM_NEGATIVE_CLASS_EXCEPTIONS
    )
    return Gram.POSITIVE if is_gram_positive else Gram.NEGATIVE


def is_gram_positive(taxonomy: Union[str, Taxonomy]) -> bool:
    """
    Determine if a bacteria is gram positive.

    Uses `bactax.gram_stain` to determine the gram stain based
    on the input information.

    Returns
    -------
    bool
        `True` if the bacteria is gram positive and `False`
        otherwise. Note that if the gram stain is `None` (because
        not enough information was given to accurately determine
        the gram stain) then `False` is returned.
    """

    return gram_stain(taxonomy) == Gram.POSITIVE


def is_gram_negative(taxonomy: Union[str, Taxonomy]) -> bool:
    """
    Determine if a bacteria is gram negative.

    Uses `bactax.gram_stain` to determine the gram stain based
    on the input information.

    Returns
    -------
    bool
        `True` if the bacteria is gram negative and `False`
        otherwise. Note that if the gram stain is `None` (because
        not enough information was given to accurately determine
        the gram stain) then `False` is returned.
    """

    return gram_stain(taxonomy) == Gram.NEGATIVE
