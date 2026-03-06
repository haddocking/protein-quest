"""Models of 3D Beacons API HUB.

Generated with
```shell
wget https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/openapi.json
uvx --from=datamodel-code-generator datamodel-codegen --input ./openapi.json \
    --output model.py  --output-model-type dataclasses.dataclass
```

After generation changes:

* ruff format
* ruff check --fix
* ruff ignores
* Trim imports
"""

# TODO document each model starting from UniprotSummary, so it shows up in API docs.
# ruff: noqa: N815
from dataclasses import dataclass
from enum import Enum


@dataclass
class AccessionListRequest:
    accessions: list[str]
    provider: str | None = None
    exclude_provider: str | None = None


class ChecksumType(Enum):
    CRC64 = "CRC64"
    MD5 = "MD5"


@dataclass
class EnsemblTranscript:
    transcript_id: str
    seqRegionStart: int
    seqRegionEnd: int
    alignment_difference: float | None = None


class EnsembleSampleFormat(Enum):
    PDB = "PDB"
    MMCIF = "MMCIF"
    BCIF = "BCIF"


@dataclass
class Entry:
    sequence: str
    checksum: str
    checksum_type: ChecksumType


class Evidence(Enum):
    EXPERIMENTALLY_DETERMINED = "EXPERIMENTALLY DETERMINED"
    COMPUTATIONAL_PREDICTED = "COMPUTATIONAL/PREDICTED"
    FROM_LITERATURE = "FROM LITERATURE"


class ExperimentalMethod1(Enum):
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    ELECTRON_MICROSCOPY = "ELECTRON MICROSCOPY"
    EPR = "EPR"
    FIBER_DIFFRACTION = "FIBER DIFFRACTION"
    FLUORESCENCE_TRANSFER = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY = "INFRARED SPECTROSCOPY"
    NEUTRON_DIFFRACTION = "NEUTRON DIFFRACTION"
    POWDER_DIFFRACTION = "POWDER DIFFRACTION"
    SOLID_STATE_NMR = "SOLID-STATE NMR"
    SOLUTION_NMR = "SOLUTION NMR"
    SOLUTION_SCATTERING = "SOLUTION SCATTERING"
    THEORETICAL_MODEL = "THEORETICAL MODEL"
    X_RAY_DIFFRACTION = "X-RAY DIFFRACTION"
    HYBRID = "HYBRID"


class FeatureType(Enum):
    CARBOHYD = "CARBOHYD"
    DOMAIN = "DOMAIN"
    ACT_SITE = "ACT_SITE"
    METAL = "METAL"
    BINDING = "BINDING"
    NON_STD = "NON_STD"
    MOD_RES = "MOD_RES"
    DISULFID = "DISULFID"
    MUTAGEN = "MUTAGEN"
    HELIX = "HELIX"
    STRAND = "STRAND"
    DISORDERED = "DISORDERED"
    INTERFACE = "INTERFACE"
    CHANNEL = "CHANNEL"


@dataclass
class HSPS:
    hsp_score: float
    hsp_bit_score: float
    hsp_align_len: int
    hsp_identity: float
    hsp_positive: float
    hsp_qseq: str
    hsp_hseq: str
    hsp_mseq: str
    hsp_expect: float


class HealthStatusEnum(Enum):
    pass_ = "pass"  # noqa: S105
    warn = "warn"
    fail = "fail"


@dataclass
class JobSubmissionErrorMessage:
    message: str | None = "Error in submitting the job, please retry!"


@dataclass
class NoJobFoundMessage:
    message: str | None = "No search job found for the given sequence, please submit the job again!"


@dataclass
class Region:
    start: int
    end: int
    annotation_value: str | None = None
    unit: str | None = None


@dataclass
class Residue:
    model_residue_label: int
    uniprot_residue_number: int
    confidence: float | None = None


@dataclass
class SearchInProgressMessage:
    message: str | None = "Search in progress, please try after sometime!"


@dataclass
class SearchSuccessMessage:
    job_id: str


@dataclass
class Seqres:
    aligned_sequence: str
    from_: int
    to: int


@dataclass
class Sequence:
    sequence: str


class SequenceIdType(Enum):
    sequence = "sequence"
    crc64 = "crc64"
    md5 = "md5"


@dataclass
class Template:
    template_id: str
    chain_id: str
    template_sequence_identity: float
    last_updated: str
    provider: str
    experimental_method: ExperimentalMethod1
    resolution: float
    preferred_assembly_id: str | None = None


@dataclass
class Uniprot:
    aligned_sequence: str
    from_: int
    to: int


@dataclass
class UniprotEntry:
    ac: str
    id: str | None = None
    uniprot_checksum: str | None = None
    sequence_length: int | None = None
    segment_start: int | None = None
    segment_end: int | None = None
    description: str | None = None


@dataclass
class ValidationError:
    loc: list[str | int]
    msg: str
    type: str


class AppSequenceSchemaConfidenceType(Enum):
    pLDDT = "pLDDT"
    QMEANDisCo = "QMEANDisCo"
    ipTM_pTM = "ipTM+pTM"


class AppSequenceSchemaEntityPolyType(Enum):
    CYCLIC_PSEUDO_PEPTIDE = "CYCLIC-PSEUDO-PEPTIDE"
    PEPTIDE_NUCLEIC_ACID = "PEPTIDE NUCLEIC ACID"
    POLYDEOXYRIBONUCLEOTIDE = "POLYDEOXYRIBONUCLEOTIDE"
    POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID = "POLYDEOXYRIBONUCLEOTIDE/POLYRIBONUCLEOTIDE HYBRID"
    POLYPEPTIDE_D_ = "POLYPEPTIDE(D)"
    POLYPEPTIDE_L_ = "POLYPEPTIDE(L)"
    POLYRIBONUCLEOTIDE = "POLYRIBONUCLEOTIDE"
    OTHER = "OTHER"


class AppSequenceSchemaEntityType(Enum):
    BRANCHED = "BRANCHED"
    MACROLIDE = "MACROLIDE"
    NON_POLYMER = "NON-POLYMER"
    POLYMER = "POLYMER"
    WATER = "WATER"


class AppSequenceSchemaExperimentalMethod(Enum):
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    ELECTRON_MICROSCOPY = "ELECTRON MICROSCOPY"
    EPR = "EPR"
    FIBER_DIFFRACTION = "FIBER DIFFRACTION"
    FLUORESCENCE_TRANSFER = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY = "INFRARED SPECTROSCOPY"
    NEUTRON_DIFFRACTION = "NEUTRON DIFFRACTION"
    X_RAY_POWDER_DIFFRACTION = "X-RAY POWDER DIFFRACTION"
    SOLID_STATE_NMR = "SOLID-STATE NMR"
    SOLUTION_NMR = "SOLUTION NMR"
    X_RAY_SOLUTION_SCATTERING = "X-RAY SOLUTION SCATTERING"
    THEORETICAL_MODEL = "THEORETICAL MODEL"
    X_RAY_DIFFRACTION = "X-RAY DIFFRACTION"
    HYBRID = "HYBRID"


class AppSequenceSchemaIdentifierCategory(Enum):
    UNIPROT = "UNIPROT"
    RFAM = "RFAM"
    CCD = "CCD"
    SMILES = "SMILES"
    INCHI = "INCHI"
    INCHIKEY = "INCHIKEY"


class AppSequenceSchemaModelCategory(Enum):
    EXPERIMENTALLY_DETERMINED = "EXPERIMENTALLY DETERMINED"
    TEMPLATE_BASED = "TEMPLATE-BASED"
    AB_INITIO = "AB-INITIO"
    CONFORMATIONAL_ENSEMBLE = "CONFORMATIONAL ENSEMBLE"


class AppSequenceSchemaModelFormat(Enum):
    PDB = "PDB"
    MMCIF = "MMCIF"
    BCIF = "BCIF"


class AppSequenceSchemaModelType(Enum):
    ATOMIC = "ATOMIC"
    DUMMY = "DUMMY"
    MIX = "MIX"


class AppSequenceSchemaOligomericState(Enum):
    MONOMER = "MONOMER"
    HOMODIMER = "HOMODIMER"
    HETERODIMER = "HETERODIMER"
    HOMO_OLIGOMER = "HOMO-OLIGOMER"
    HETERO_OLIGOMER = "HETERO-OLIGOMER"


class AppUniprotSchemaConfidenceType(Enum):
    pLDDT = "pLDDT"
    QMEANDisCo = "QMEANDisCo"
    ipTM_pTM = "ipTM+pTM"


class AppUniprotSchemaEntityPolyType(Enum):
    CYCLIC_PSEUDO_PEPTIDE = "CYCLIC-PSEUDO-PEPTIDE"
    PEPTIDE_NUCLEIC_ACID = "PEPTIDE NUCLEIC ACID"
    POLYDEOXYRIBONUCLEOTIDE = "POLYDEOXYRIBONUCLEOTIDE"
    POLYDEOXYRIBONUCLEOTIDE_POLYRIBONUCLEOTIDE_HYBRID = "POLYDEOXYRIBONUCLEOTIDE/POLYRIBONUCLEOTIDE HYBRID"
    POLYPEPTIDE_D_ = "POLYPEPTIDE(D)"
    POLYPEPTIDE_L_ = "POLYPEPTIDE(L)"
    POLYRIBONUCLEOTIDE = "POLYRIBONUCLEOTIDE"
    OTHER = "OTHER"


class AppUniprotSchemaEntityType(Enum):
    BRANCHED = "BRANCHED"
    MACROLIDE = "MACROLIDE"
    NON_POLYMER = "NON-POLYMER"
    POLYMER = "POLYMER"
    WATER = "WATER"


class AppUniprotSchemaExperimentalMethod(Enum):
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    ELECTRON_MICROSCOPY = "ELECTRON MICROSCOPY"
    EPR = "EPR"
    FIBER_DIFFRACTION = "FIBER DIFFRACTION"
    FLUORESCENCE_TRANSFER = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY = "INFRARED SPECTROSCOPY"
    NEUTRON_DIFFRACTION = "NEUTRON DIFFRACTION"
    X_RAY_POWDER_DIFFRACTION = "X-RAY POWDER DIFFRACTION"
    SOLID_STATE_NMR = "SOLID-STATE NMR"
    SOLUTION_NMR = "SOLUTION NMR"
    X_RAY_SOLUTION_SCATTERING = "X-RAY SOLUTION SCATTERING"
    THEORETICAL_MODEL = "THEORETICAL MODEL"
    X_RAY_DIFFRACTION = "X-RAY DIFFRACTION"
    HYBRID = "HYBRID"


class AppUniprotSchemaIdentifierCategory(Enum):
    UNIPROT = "UNIPROT"
    RFAM = "RFAM"
    CCD = "CCD"
    SMILES = "SMILES"
    INCHI = "INCHI"
    INCHIKEY = "INCHIKEY"


class AppUniprotSchemaModelCategory(Enum):
    EXPERIMENTALLY_DETERMINED = "EXPERIMENTALLY DETERMINED"
    TEMPLATE_BASED = "TEMPLATE-BASED"
    AB_INITIO = "AB-INITIO"
    CONFORMATIONAL_ENSEMBLE = "CONFORMATIONAL ENSEMBLE"


class AppUniprotSchemaModelFormat(Enum):
    PDB = "PDB"
    MMCIF = "MMCIF"
    BCIF = "BCIF"


class AppUniprotSchemaModelType(Enum):
    ATOMIC = "ATOMIC"
    DUMMY = "DUMMY"
    MIX = "MIX"


class AppUniprotSchemaOligomericState(Enum):
    MONOMER = "MONOMER"
    HOMODIMER = "HOMODIMER"
    HETERODIMER = "HETERODIMER"
    HOMO_OLIGOMER = "HOMO-OLIGOMER"
    HETERO_OLIGOMER = "HETERO-OLIGOMER"


@dataclass
class FeatureItem:
    type: FeatureType
    description: str
    source_name: str | None = None
    source_url: str | None = None
    evidence: Evidence | None = None
    residues: list[int] | None = None
    regions: list[Region] | None = None


@dataclass
class HTTPValidationError:
    detail: list[ValidationError] | None = None


@dataclass
class HealthStatus:
    status: HealthStatusEnum
    service_id: str
    beacons_api_version: str
    output: str | None = None


@dataclass
class Segment:
    seqres: Seqres
    uniprot: Uniprot
    residues: list[Residue]
    templates: list[Template] | None = None


@dataclass
class AppSequenceSchemaEntity:
    entity_type: AppSequenceSchemaEntityType
    description: str
    chain_ids: list[str]
    entity_poly_type: AppSequenceSchemaEntityPolyType | None = None
    identifier: str | None = None
    identifier_category: AppSequenceSchemaIdentifierCategory | None = None


@dataclass
class AppSequenceSchemaSummaryItems:
    model_identifier: str
    model_category: AppSequenceSchemaModelCategory
    model_url: str
    model_format: AppSequenceSchemaModelFormat
    provider: str
    created: str
    sequence_identity: float
    coverage: float
    entities: list[AppSequenceSchemaEntity]
    model_type: AppSequenceSchemaModelType | None = None
    model_page_url: str | None = None
    number_of_conformers: float | None = None
    ensemble_sample_url: str | None = None
    ensemble_sample_format: AppSequenceSchemaModelFormat | None = None
    experimental_method: AppSequenceSchemaExperimentalMethod | None = None
    resolution: float | None = None
    confidence_type: AppSequenceSchemaConfidenceType | None = None
    confidence_version: str | None = None
    confidence_avg_local_score: float | None = None
    oligomeric_state: AppSequenceSchemaOligomericState | None = None
    oligomeric_state_confidence: float | None = None
    preferred_assembly_id: str | None = None


@dataclass
class AppUniprotSchemaEntity:
    entity_type: AppUniprotSchemaEntityType
    description: str
    chain_ids: list[str]
    entity_poly_type: AppUniprotSchemaEntityPolyType | None = None
    identifier: str | None = None
    identifier_category: AppUniprotSchemaIdentifierCategory | None = None


@dataclass
class AppUniprotSchemaSummaryItems:
    model_identifier: str
    model_category: AppUniprotSchemaModelCategory
    model_url: str
    model_format: AppUniprotSchemaModelFormat
    provider: str
    created: str
    sequence_identity: float
    uniprot_start: int
    uniprot_end: int
    coverage: float
    entities: list[AppUniprotSchemaEntity]
    model_type: AppUniprotSchemaModelType | None = None
    model_page_url: str | None = None
    number_of_conformers: float | None = None
    ensemble_sample_url: str | None = None
    ensemble_sample_format: EnsembleSampleFormat | None = None
    experimental_method: AppUniprotSchemaExperimentalMethod | None = None
    resolution: float | None = None
    confidence_type: AppUniprotSchemaConfidenceType | None = None
    confidence_version: str | None = None
    confidence_avg_local_score: float | None = None
    oligomeric_state: AppUniprotSchemaOligomericState | None = None
    oligomeric_state_confidence: float | None = None
    preferred_assembly_id: str | None = None


@dataclass
class Annotation:
    accession: str
    sequence: str
    id: str | None = None
    annotation: list[FeatureItem] | None = None


@dataclass
class Chain:
    chain_id: str
    segments: list[Segment] | None = None


type Chains = list[Chain]


@dataclass
class Detailed:
    summary: AppUniprotSchemaSummaryItems
    chains: Chains


type HealthResponse = list[HealthStatus]


@dataclass
class Overview:
    summary: AppUniprotSchemaSummaryItems


@dataclass
class SequenceOverview:
    summary: AppSequenceSchemaSummaryItems


@dataclass
class SequenceSummary:
    entry: Entry
    structures: list[SequenceOverview]


@dataclass
class UniprotDetails:
    uniprot_entry: UniprotEntry | None = None
    structures: list[Detailed] | None = None


@dataclass
class UniprotSummary:
    uniprot_entry: UniprotEntry | None = None
    structures: list[Overview] | None = None


@dataclass
class SearchAccession:
    accession: str
    id: str
    description: str
    hit_length: int
    hit_hsps: list[HSPS]
    hit_uni_ox: int
    hit_uni_os: str
    hit_com_os: str
    title: str
    summary: UniprotSummary | None = None


@dataclass
class UniprotMapping:
    ensembl_transcript: EnsemblTranscript
    uniprot_accession: UniprotSummary


@dataclass
class EnsemblSummary:
    ensembl_id: str
    species: str
    taxid: str
    uniprot_mappings: list[UniprotMapping]
