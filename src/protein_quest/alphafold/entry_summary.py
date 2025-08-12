# ruff: noqa: N815 allow camelCase follow what api returns
from dataclasses import dataclass


@dataclass
class EntrySummary:
    """Dataclass representing a summary of an AlphaFold entry.

    Modelled after EntrySummary in https://alphafold.ebi.ac.uk/api/openapi.json
    """

    entryId: str
    uniprotAccession: str
    uniprotId: str
    uniprotDescription: str
    taxId: int
    organismScientificName: str
    uniprotStart: int
    uniprotEnd: int
    uniprotSequence: str
    modelCreatedDate: str
    latestVersion: int
    allVersions: list[int]
    bcifUrl: str
    cifUrl: str
    pdbUrl: str
    paeImageUrl: str
    paeDocUrl: str
    gene: str | None = None
    sequenceChecksum: str | None = None
    sequenceVersionDate: str | None = None
    amAnnotationsUrl: str | None = None
    amAnnotationsHg19Url: str | None = None
    amAnnotationsHg38Url: str | None = None
    isReviewed: bool | None = None
    isReferenceProteome: bool | None = None
