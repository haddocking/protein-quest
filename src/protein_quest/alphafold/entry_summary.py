# ruff: noqa: N815 allow camelCase follow what api returns
from dataclasses import dataclass

from yarl import URL


@dataclass
class EntrySummary:
    """Dataclass representing a summary of an AlphaFold entry.

    Modelled after NewEntrySummary in [https://alphafold.ebi.ac.uk/api/openapi.json](https://alphafold.ebi.ac.uk/api/openapi.json)
    with URL types and without deprecated fields.
    """

    allVersions: list[int]
    bcifUrl: URL
    cifUrl: URL
    entityType: str
    fractionPlddtConfident: float
    fractionPlddtLow: float
    fractionPlddtVeryHigh: float
    fractionPlddtVeryLow: float
    globalMetricValue: float
    isUniProt: bool
    latestVersion: int
    modelCreatedDate: str
    modelEntityId: str
    paeDocUrl: URL
    pdbUrl: URL
    providerId: str
    sequence: str
    sequenceChecksum: str
    sequenceEnd: int
    sequenceStart: int
    sequenceVersionDate: str
    toolUsed: str
    alternativeNames: list[str] | None = None
    amAnnotationsHg19Url: URL | None = None
    amAnnotationsHg38Url: URL | None = None
    amAnnotationsUrl: URL | None = None
    catalyticActivities: list[str] | None = None
    complexName: str | None = None
    functions: list[str] | None = None
    gene: str | None = None
    geneSynonyms: list[str] | None = None
    ipSAE: float | None = None
    ipTM: float | None = None
    isUniProtReferenceProteome: bool | None = None
    isUniProtReviewed: bool | None = None
    keywords: list[str] | None = None
    msaUrl: URL | None = None
    organismCommonNames: list[str] | None = None
    organismScientificName: str | None = None
    organismSynonyms: list[str] | None = None
    plddtDocUrl: URL | None = None
    proteinFullNames: list[str] | None = None
    proteinShortNames: list[str] | None = None
    stoichiometry: int | None = None
    taxId: int | None = None
    taxonomyLineage: list[str] | None = None
    # uniprotAccession is isoform id (<uniprot_accession>-<isoform number>) when entry has multiple isoforms.
    uniprotAccession: str | None = None
    uniprotDescription: str | None = None
    uniprotId: str | None = None
