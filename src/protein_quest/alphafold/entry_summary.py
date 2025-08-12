# ruff: noqa: N815 allow camelCase follow what api returns
from dataclasses import dataclass


@dataclass
class EntrySummary:
    """Dataclass representing a summary of an AlphaFold entry.

    Modelled after EntrySummary in https://alphafold.ebi.ac.uk/api/openapi.json
    """

    entryId: str
    gene: str | None
    sequenceChecksum: str | None
    sequenceVersionDate: str | None
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
    amAnnotationsUrl: str | None
    amAnnotationsHg19Url: str | None
    amAnnotationsHg38Url: str | None
    isReviewed: bool | None
    isReferenceProteome: bool | None
