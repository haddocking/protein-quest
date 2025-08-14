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
    # TODO add new fields from https://alphafold.ebi.ac.uk/#/public-api/get_uniprot_summary_api_uniprot_summary__qualifier__json_get
    # TODO like fractionPlddt* fields which can be used in filter_files_on_confidence()
