from typing import Union, Iterable, Optional, Sequence
import os
from pathlib import Path
from zipfile import ZipFile
import fnmatch
import hashlib
import warnings

from pydantic import BaseModel
from strenum import StrEnum

from bates_reconstruction_data.utils import download


NCBI_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"


class AnnotationType(StrEnum):
    GENOME_GFF = "GENOME_GFF"
    GENOME_GBFF = "GENOME_GBFF"
    RNA_FASTA = "RNA_FASTA"
    PROT_FASTA = "PROT_FASTA"
    GENOME_GTF = "GENOME_GTF"
    CDS_FASTA = "CDS_FASTA"
    GENOME_FASTA = "GENOME_FASTA"
    SEQUENCE_REPORT = "SEQUENCE_REPORT"


class Hydration(StrEnum):
    FULLY_HYDRATED = "FULLY_HYDRATED"
    DATA_REPORT_ONLY = "DATA_REPORT_ONLY"


class GenomeDownloadRequest(BaseModel):

    accessions: Sequence[str]
    chromosomes: Optional[Sequence[str]] = None
    include_annotation_type: Optional[Sequence[AnnotationType]] = None
    hydrated: Optional[Hydration] = None
    include_tsv: bool = False


class NCBI:

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.genome = GenomeClient()
        self.gene = GeneClient()


class GenomeClient:

    def __init__(self):
        ...

    def download(
            self,
            accessions: Union[str, Iterable[str]],
            path: Union[str, os.PathLike],
            *,
            chromosomes: Optional[Union[str, Iterable[str]]] = None,
            include_annotation_type: Optional[Union[AnnotationType, Iterable[AnnotationType]]] = None,
            hydrated: Hydration = Hydration.FULLY_HYDRATED,
            include_tsv: bool = False,
            validate: bool = False
        ) -> "NCBIDataPackage":
        """

        """

        if isinstance(accessions, str):
            accessions = [accessions]
        
        if isinstance(include_annotation_type, str): # True for both str and AnnotationType (which is a StrEnum)
            include_annotation_type = [AnnotationType(include_annotation_type)]

        model = GenomeDownloadRequest(
            accessions=accessions,
            chromosomes=chromosomes,
            include_annotation_type=include_annotation_type,
            hydrated=hydrated,
            include_tsv=include_tsv
        )

        url = f"{NCBI_BASE_URL}/genome/download"
        download.download(
            url,
            path,
            method="POST",
            json=model.model_dump(),
            headers={"content-type": "application/json"}
        )

        return NCBIDataPackage(path, validate=validate)


class GeneClient:
    pass


class NCBIDataPackage:
    """
    Current idea: make this kind of an opinionated wrapper around a ZipFile that makes
    pulling data out of an NCBI dataset archive easier.
    """
    
    def __init__(self, path: os.PathLike, *, scan: bool = True, validate: bool = False):
        self.path = Path(path)
        self.datasets: list[NCBIDataset] = []
        
        if scan:
            self.scan()

        if validate and not self.is_valid():
            warnings.warn(f"{self.path} failed the md5 checksum validation.")

    def scan(self):
        fn_pattern = os.path.join("ncbi_dataset", "data", "*", "*")
        accessions = set()
        with ZipFile(self.path) as z:
            for filename in z.namelist():
                if fnmatch.fnmatch(filename, fn_pattern):
                    accession = filename.split(os.path.sep)[-2]
                    if accession not in accessions:
                        accessions.add(accession)
                        self.datasets.append(NCBIDataset(self.path, accession))

    def is_valid(self) -> bool:
        """
        Checks the package contents against the md5 checksum that is included in
        the package.
        """
        with ZipFile(self.path) as z:
            checksums = z.read("md5sum.txt").decode("utf-8")
            for line in checksums.splitlines():
                checksum, filename = line.split()
                with z.open(filename) as file:
                    hash = hashlib.md5()
                    for chunk in file:
                        hash.update(chunk)
                    if hash.hexdigest() != checksum:
                        return False
        return True
    
    def zipfile(self) -> ZipFile:
        return ZipFile(self.path)


class NCBIDataset:

    def __init__(self, archive_path: os.PathLike, accession: str):
        self.data_package = NCBIDataPackage(archive_path, scan=False)
        self.accession = accession

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({repr(self.data_package.path)}, {repr(self.accession)})"
    
    def list_files(self) -> list[str]:
        files = []
        fn_pattern = os.path.join("ncbi_dataset", "data", self.accession, "*")
        with ZipFile(self.data_package.path) as z:
            for file in z.namelist():
                if fnmatch.fnmatch(file, fn_pattern):
                    filename = file.split(os.path.sep)[-1]
                    files.append(os.path.join(self.accession, filename))
        return files
