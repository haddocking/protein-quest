"""Models of 3D Beacons API HUB.

Generated with
```shell
wget https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/v2/openapi.json
uvx --from=datamodel-code-generator datamodel-codegen \
    --output-model-type dataclasses.dataclass \
    --enum-field-as-literal all --disable-future-imports \
    --use-annotated --use-field-description \
    --openapi-scopes schemas parameters paths --include-path-parameters \
    --input ./openapi.json --output model.py
```
When 3D beaconds API changes this module should be regenerated and the above command and steps below should be repeated.

After generation changes:

* ruff format
* ruff check --fix
* ruff ignores
* extract Provider type

Also document UniprotSummary and its children:

* Use google style docstrings using `Attributes:` section
* as mkdocs does not understand docstring per property
* make llm do it with following prompts:
```
Document the children of UniprotSummary and AccessionListRequest.
If class has attributes with their own docstrings move them to
the `Attributes:` section of the classes docstring.
Do not add types in the docstring as they are already defined in the code.
```
and
```
Add docstring to types that do not have any.
Convert name of type from CamelCase to Camel case as docstring.
```

"""

# ruff: noqa: N815
from dataclasses import dataclass
from typing import Any, Literal


@dataclass
class AccessionListRequest:
    """Accession list request

    Attributes:
        accessions: A list of UniProt accessions
        provider: Name of the model provider
        exclude_provider: Provider to exclude.

    """

    accessions: list[str]
    provider: str | None = None
    exclude_provider: str | None = None


type ChecksumType = Literal["CRC64", "MD5"]
"""
Type of checksum enumeration
"""


type AlignmentDifference = float
"""
Fraction of alignment difference (0 to 1)
"""


@dataclass
class EnsemblTranscript:
    transcript_id: str
    """
    Transcript identifier
    """
    seqRegionStart: int
    """
    Start position of the transcript
    """
    seqRegionEnd: int
    """
    End position of the transcript
    """
    alignment_difference: AlignmentDifference | None = None
    """
    Fraction of alignment difference (0 to 1)
    """


type EnsembleSampleFormat = Literal["PDB", "MMCIF", "BCIF"]
"""
Ensemble Sample Format
"""


@dataclass
class Entry:
    sequence: str
    """
    The protein sequence
    """
    checksum: str
    """
    CRC64 or MD5 checksum of the sequence
    """
    checksum_type: ChecksumType
    """
    Type of the checksum
    """


type Evidence = Literal["EXPERIMENTALLY DETERMINED", "COMPUTATIONAL/PREDICTED", "FROM LITERATURE"]
"""
Evidence
"""


type ExperimentalMethod1 = Literal[
    "ELECTRON CRYSTALLOGRAPHY",
    "ELECTRON MICROSCOPY",
    "EPR",
    "FIBER DIFFRACTION",
    "FLUORESCENCE TRANSFER",
    "INFRARED SPECTROSCOPY",
    "NEUTRON DIFFRACTION",
    "POWDER DIFFRACTION",
    "SOLID-STATE NMR",
    "SOLUTION NMR",
    "SOLUTION SCATTERING",
    "THEORETICAL MODEL",
    "X-RAY DIFFRACTION",
    "HYBRID",
]
"""
Experimental Method 1
"""


type FeatureType = Literal[
    "CARBOHYD",
    "DOMAIN",
    "ACT_SITE",
    "METAL",
    "BINDING",
    "NON_STD",
    "MOD_RES",
    "DISULFID",
    "MUTAGEN",
    "HELIX",
    "STRAND",
    "DISORDERED",
    "INTERFACE",
    "CHANNEL",
]
"""
Feature Type
"""


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


type HealthStatusEnum = Literal["pass", "warn", "fail"]
"""
Health status enumeration based on the OpenAPI spec
"""


@dataclass
class JobSubmissionErrorMessage:
    message: str | None = "Error in submitting the job, please retry!"


@dataclass
class NoJobFoundMessage:
    message: str | None = "No search job found for the given sequence, please submit the job again!"


@dataclass
class Region:
    start: int
    """
    The first position of the annotation
    """
    end: int
    """
    The last position of the annotation
    """
    annotation_value: str | None = None
    """
    The value of the annotation
    """
    unit: str | None = None
    """
    The unit of the annotation value, if applicable
    """


@dataclass
class Residue:
    model_residue_label: int
    """
    Model residue index
    """
    uniprot_residue_number: int
    """
    UniProt residue index
    """
    confidence: float | None = None
    """
    Confidence score in the range of [0,1]
    """


@dataclass
class SearchInProgressMessage:
    message: str | None = "Search in progress, please try after sometime!"


@dataclass
class SearchSuccessMessage:
    job_id: str


@dataclass
class Seqres:
    aligned_sequence: str
    """
    Sequence of the model
    """
    from_: int
    """
    1-indexed first residue
    """
    to: int
    """
    1-indexed last residue
    """


@dataclass
class Sequence:
    sequence: str


type SequenceIdType = Literal["sequence", "crc64", "md5"]
"""
Sequence Id Type
"""


@dataclass
class Template:
    template_id: str
    """
    Identifier of the template
    """
    chain_id: str
    """
    Identifier of the chain of the template; this is label_asym_id in mmCIF
    """
    template_sequence_identity: float
    """
    Sequence identity of the template with the  UniProt accession, in the range of [0,1]

    """
    last_updated: str
    """
    Date of release of the last update in  the format of YYYY-MM-DD

    """
    provider: str
    """
    Provider of the template
    """
    experimental_method: ExperimentalMethod1
    """
    Experimental method used to determine the template
    """
    resolution: float
    """
    Resolution of the template, in Angstrom
    """
    preferred_assembly_id: str | None = None
    """
    Identifier of the preferred assembly of the template
    """


@dataclass
class Uniprot:
    aligned_sequence: str
    """
    Sequence of the UniProt accession
    """
    from_: int
    """
    1-indexed first residue
    """
    to: int
    """
    1-indexed last residue
    """


@dataclass
class UniprotEntry:
    """UniProt entry

    Attributes:
        ac: UniProt accession.
        id: UniProt identifier.
        uniprot_checksum: CRC64 checksum of the UniProt sequence.
        sequence_length: Length of the UniProt sequence.
        segment_start: 1-indexed first residue of the UniProt sequence segment.
        segment_end: 1-indexed last residue of the UniProt sequence segment.
        description: Description of the UniProt entry.
    """

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


type AppSequenceSchemaConfidenceType = Literal["pLDDT", "QMEANDisCo", "ipTM+pTM"]
"""
Confidence type enumeration
"""


type AppSequenceSchemaEntityPolyType = Literal[
    "CYCLIC-PSEUDO-PEPTIDE",
    "PEPTIDE NUCLEIC ACID",
    "POLYDEOXYRIBONUCLEOTIDE",
    "POLYDEOXYRIBONUCLEOTIDE/POLYRIBONUCLEOTIDE HYBRID",
    "POLYPEPTIDE(D)",
    "POLYPEPTIDE(L)",
    "POLYRIBONUCLEOTIDE",
    "OTHER",
]
"""
Entity poly type enumeration
"""


type AppSequenceSchemaEntityType = Literal["BRANCHED", "MACROLIDE", "NON-POLYMER", "POLYMER", "WATER"]
"""
Entity type enumeration
"""


type AppSequenceSchemaExperimentalMethod = Literal[
    "ELECTRON CRYSTALLOGRAPHY",
    "ELECTRON MICROSCOPY",
    "EPR",
    "FIBER DIFFRACTION",
    "FLUORESCENCE TRANSFER",
    "INFRARED SPECTROSCOPY",
    "NEUTRON DIFFRACTION",
    "X-RAY POWDER DIFFRACTION",
    "SOLID-STATE NMR",
    "SOLUTION NMR",
    "X-RAY SOLUTION SCATTERING",
    "THEORETICAL MODEL",
    "X-RAY DIFFRACTION",
    "HYBRID",
]
"""
Experimental method enumeration
"""


type AppSequenceSchemaIdentifierCategory = Literal["UNIPROT", "RFAM", "CCD", "SMILES", "INCHI", "INCHIKEY"]
"""
Identifier category enumeration
"""


type AppSequenceSchemaModelCategory = Literal[
    "EXPERIMENTALLY DETERMINED",
    "TEMPLATE-BASED",
    "AB-INITIO",
    "CONFORMATIONAL ENSEMBLE",
]
"""
Model category enumeration
"""


type AppSequenceSchemaModelFormat = Literal["PDB", "MMCIF", "BCIF"]
"""
Model format enumeration
"""


type AppSequenceSchemaModelType = Literal["ATOMIC", "DUMMY", "MIX"]
"""
Model type enumeration
"""


type AppSequenceSchemaOligomericState = Literal[
    "MONOMER", "HOMODIMER", "HETERODIMER", "HOMO-OLIGOMER", "HETERO-OLIGOMER"
]
"""
Oligomeric state enumeration
"""


type Resolution = float
"""
Resolution of the model in Angstrom
"""


type AppUniprotSchemaConfidenceType = Literal["pLDDT", "QMEANDisCo", "ipTM+pTM"]
"""
App Uniprot Schema Confidence Type
"""


type AppUniprotSchemaEntityPolyType = Literal[
    "CYCLIC-PSEUDO-PEPTIDE",
    "PEPTIDE NUCLEIC ACID",
    "POLYDEOXYRIBONUCLEOTIDE",
    "POLYDEOXYRIBONUCLEOTIDE/POLYRIBONUCLEOTIDE HYBRID",
    "POLYPEPTIDE(D)",
    "POLYPEPTIDE(L)",
    "POLYRIBONUCLEOTIDE",
    "OTHER",
]
"""
App Uniprot Schema Entity Poly Type
"""


type AppUniprotSchemaEntityType = Literal["BRANCHED", "MACROLIDE", "NON-POLYMER", "POLYMER", "WATER"]
"""
App Uniprot Schema Entity Type
"""


type AppUniprotSchemaExperimentalMethod = Literal[
    "ELECTRON CRYSTALLOGRAPHY",
    "ELECTRON MICROSCOPY",
    "EPR",
    "FIBER DIFFRACTION",
    "FLUORESCENCE TRANSFER",
    "INFRARED SPECTROSCOPY",
    "NEUTRON DIFFRACTION",
    "X-RAY POWDER DIFFRACTION",
    "SOLID-STATE NMR",
    "SOLUTION NMR",
    "X-RAY SOLUTION SCATTERING",
    "THEORETICAL MODEL",
    "X-RAY DIFFRACTION",
    "HYBRID",
]
"""
App Uniprot Schema Experimental Method
"""


type AppUniprotSchemaIdentifierCategory = Literal["UNIPROT", "RFAM", "CCD", "SMILES", "INCHI", "INCHIKEY"]
"""
App Uniprot Schema Identifier Category
"""


type AppUniprotSchemaModelCategory = Literal[
    "EXPERIMENTALLY DETERMINED",
    "TEMPLATE-BASED",
    "AB-INITIO",
    "CONFORMATIONAL ENSEMBLE",
]
"""
App Uniprot Schema Model Category
"""


type AppUniprotSchemaModelFormat = Literal["PDB", "MMCIF", "BCIF"]
"""
App Uniprot Schema Model Format
"""


type AppUniprotSchemaModelType = Literal["ATOMIC", "DUMMY", "MIX"]
"""
App Uniprot Schema Model Type
"""


type AppUniprotSchemaOligomericState = Literal[
    "MONOMER", "HOMODIMER", "HETERODIMER", "HOMO-OLIGOMER", "HETERO-OLIGOMER"
]
"""
App Uniprot Schema Oligomeric State
"""

Provider = Literal[
    "pdbe",
    "ped",
    "swissmodel",
    "alphafold",
    "sasbdb",
    "alphafill",
    "hegelab",
    "modelarchive",
    "isoformio",
    "levylab",
]
"""
Provider
"""


@dataclass
class UniprotSummaryQualifierJsonGetParameters:
    qualifier: Any
    """
    UniProtKB accession number (AC), entry name (ID) or CRC64 checksum of the UniProt sequence
    """
    provider: Provider | None = None
    template: Any = None
    """
    Template is 4 letter PDB code, or 4 letter code with assembly ID and chain for SMTL entries
    """
    range: Any = None
    """
    Specify a UniProt sequence residue range
    """
    exclude_provider: Provider | None = None
    """
    Provider to exclude.
    """
    uniprot_checksum: str | None = None
    """
    CRC64 checksum of the UniProt sequence
    """


@dataclass
class UniprotQualifierJsonGetParameters:
    qualifier: Any
    """
    UniProtKB accession number (AC), entry name (ID) or CRC64 checksum of the UniProt sequence
    """
    provider: Literal["swissmodel", "pdbe"] | None = None
    template: Any = None
    """
    Template is 4 letter PDB code, or 4 letter code with assembly ID and chain for SMTL entries
    """
    range: Any = None
    """
    Specify a UniProt sequence residue range
    """
    uniprot_checksum: str | None = None
    """
    CRC64 checksum of the UniProt sequence
    """


@dataclass
class EnsemblSummaryQualifierJsonGetParameters:
    qualifier: Any
    """
    Ensembl identifier.
    """
    provider: Provider | None = None


@dataclass
class SequenceResultGetParameters:
    job_id: str


@dataclass
class SequenceSummaryGetParameters:
    id: str
    """
    Identifier for the type specified in the type parameter
    """
    type: SequenceIdType | None = "sequence"
    """
    Type of the identifier
    """


@dataclass
class AnnotationsUniprotQualifierJsonGetParameters:
    uniprot_qualifier: Any
    """
    UniProtKB accession number (AC) or entry name (ID)
    """
    type: FeatureType
    """
    Annotation type
    """
    provider: Literal["pdbe", "alphafold"] | None = None
    range: Any | None = None
    """
    Specify a UniProt sequence residue range; separated by hyphen(-).
    """


@dataclass
class FeatureItem:
    type: FeatureType
    """
    Type of the annotation
    """
    description: str
    """
    Description/Label of the annotation
    """
    source_name: str | None = None
    """
    Name of the source of the annotations, i.e. where is the data from
    """
    source_url: str | None = None
    """
    URL of the source of the annotation, i.e. where to find more data
    """
    evidence: Evidence | None = None
    """
    Evidence category of the annotation
    """
    residues: list[int] | None = None
    """
    An array of residue indices
    """
    regions: list[Region] | None = None


@dataclass
class HTTPValidationError:
    detail: list[ValidationError] | None = None


@dataclass
class HealthStatus:
    status: HealthStatusEnum
    """
    The status of the service
    """
    service_id: str
    """
    The identifier of the service, corresponds to the provider ID in the registry
    """
    beacons_api_version: str
    """
    The version of the 3DBeacons API
    """
    output: str | None = None
    """
    Raw error output, in case of 'fail' or 'warn' states
    """


@dataclass
class Segment:
    seqres: Seqres
    """
    Information on the sequence of the model
    """
    uniprot: Uniprot
    residues: list[Residue]
    templates: list[Template] | None = None
    """
    Information on the template(s) used for the model
    """


@dataclass
class AppSequenceSchemaEntity:
    entity_type: AppSequenceSchemaEntityType
    """
    Type of the molecular entity
    """
    description: str
    """
    Textual label of the molecule
    """
    chain_ids: list[str]
    """
    List of chain identifiers of the molecule
    """
    entity_poly_type: AppSequenceSchemaEntityPolyType | None = None
    """
    Type of the molecular entity polymer
    """
    identifier: str | None = None
    """
    Identifier of the molecule
    """
    identifier_category: AppSequenceSchemaIdentifierCategory | None = None
    """
    Category of the identifier
    """


@dataclass
class AppSequenceSchemaSummaryItems:
    model_identifier: str
    """
    Identifier of the model, such as PDB id
    """
    model_category: AppSequenceSchemaModelCategory
    """
    Category of the model
    """
    model_url: str
    """
    URL of the model coordinates
    """
    model_format: AppSequenceSchemaModelFormat
    """
    File format of the coordinates
    """
    provider: str
    """
    Name of the model provider
    """
    created: str
    """
    Date of release of model generation in the format of YYYY-MM-DD
    """
    sequence_identity: float
    """
    Sequence identity of model to UniProt sequence
    """
    coverage: float
    """
    Fraction of UniProt sequence covered by the model
    """
    entities: list[AppSequenceSchemaEntity]
    """
    List of molecular entities in the model
    """
    model_type: AppSequenceSchemaModelType | None = None
    """
    Defines if coordinates are atomic-level or contain dummy atoms
    """
    model_page_url: str | None = None
    """
    URL of a web page showing the model
    """
    number_of_conformers: float | None = None
    """
    Number of conformers in a conformational ensemble
    """
    ensemble_sample_url: str | None = None
    """
    URL of a sample of conformations from ensemble
    """
    ensemble_sample_format: AppSequenceSchemaModelFormat | None = None
    """
    File format of the sample coordinates
    """
    experimental_method: AppSequenceSchemaExperimentalMethod | None = None
    """
    Experimental method used to determine structure
    """
    resolution: Resolution | None = None
    """
    Resolution of the model in Angstrom
    """
    confidence_type: AppSequenceSchemaConfidenceType | None = None
    """
    Type of confidence measure
    """
    confidence_version: str | None = None
    """
    Version of confidence measure software
    """
    confidence_avg_local_score: float | None = None
    """
    Average of confidence measures
    """
    oligomeric_state: AppSequenceSchemaOligomericState | None = None
    """
    Oligomeric state of the model
    """
    oligomeric_state_confidence: float | None = None
    """
    Confidence in oligomeric state
    """
    preferred_assembly_id: str | None = None
    """
    Identifier of preferred assembly
    """


@dataclass
class AppUniprotSchemaEntity:
    """Molecular entity in a UniProt structure model.

    Attributes:
        entity_type: Type of the molecular entity; similar to `_entity.type` in mmCIF.
        description: Textual label of the molecule.
        chain_ids: List of label_asym identifiers (`chain_id` in PDB format) for the molecule.
        entity_poly_type: Type of the molecular entity polymer; similar to `_entity_poly.type` in mmCIF.
        identifier: Identifier of the molecule.
        identifier_category: Category of the identifier.
    """

    entity_type: AppUniprotSchemaEntityType
    description: str
    chain_ids: list[str]
    entity_poly_type: AppUniprotSchemaEntityPolyType | None = None
    identifier: str | None = None
    identifier_category: AppUniprotSchemaIdentifierCategory | None = None


@dataclass
class AppUniprotSchemaSummaryItems:
    """Summary information for a UniProt-associated structure model.

    Attributes:
        model_identifier: Identifier of the model, such as PDB ID.
        model_category: Category of the model.
        model_url: URL of the model coordinates.
        model_format: File format of the coordinates.
        provider: Name of the model provider.
        created: Date of model release in `YYYY-MM-DD` format.
        sequence_identity: Sequence identity of the model to the UniProt sequence in range `[0, 1]`.
        uniprot_start: 1-indexed first residue of the model in UniProt numbering.
        uniprot_end: 1-indexed last residue of the model in UniProt numbering.
        coverage: Fraction of UniProt sequence covered by the model in range `[0, 1]`.
        entities: List of molecular entities in the model.
        model_type: Whether coordinates are atomic, dummy atoms, or mixed.
        model_page_url: Provider web page URL for the model.
        number_of_conformers: Number of conformers in a conformational ensemble.
        ensemble_sample_url: URL of a sample of conformations from an ensemble.
        ensemble_sample_format: File format of the sample coordinates.
        experimental_method: Experimental method used to determine the structure.
        resolution: Resolution of the model in Angstrom.
        confidence_type: Type of confidence measure.
        confidence_version: Version of the confidence measure software.
        confidence_avg_local_score: Average confidence score.
        oligomeric_state: Oligomeric state of the model.
        oligomeric_state_confidence: Confidence in the oligomeric state assignment.
        preferred_assembly_id: Identifier of the preferred assembly in the model.
    """

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
    """
    A UniProt accession
    """
    sequence: str
    """
    The sequence of the protein
    """
    id: str | None = None
    """
    A UniProt identifier
    """
    annotation: list[FeatureItem] | None = None


@dataclass
class Chain:
    chain_id: str
    segments: list[Segment] | None = None


type Chains = list[Chain]
"""
Chains
"""


@dataclass
class Detailed:
    summary: AppUniprotSchemaSummaryItems
    chains: Chains


type HealthResponse = list[HealthStatus]
"""
Response model for /health endpoint - returns array of status objects
"""


@dataclass
class Overview:
    """Overview of a UniProt-associated structure.

    Attributes:
        summary: Summary information for the structure.
    """

    summary: AppUniprotSchemaSummaryItems


@dataclass
class SequenceOverview:
    summary: AppSequenceSchemaSummaryItems
    """
    Summary information for the structure
    """


@dataclass
class SequenceSummary:
    entry: Entry
    """
    Entry information including sequence and checksum
    """
    structures: list[SequenceOverview]
    """
    List of available structures
    """


@dataclass
class UniprotDetails:
    uniprot_entry: UniprotEntry | None = None
    structures: list[Detailed] | None = None


@dataclass
class UniprotSummary:
    """Uniprot summary

    Attributes:
        uniprot_entry: Uniprot entry
        structures: Structures
    """

    uniprot_entry: UniprotEntry | None = None
    structures: list[Overview] | None = None


type UniprotSummaryPostResponse = list[UniprotSummary]
"""
Uniprot Summary Post Response
"""


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


type SequenceResultGetResponse = list[SearchAccession]
"""
Sequence Result Get Response
"""


@dataclass
class EnsemblSummary:
    ensembl_id: str
    """
    An Ensembl identifier
    """
    species: str
    """
    Species name
    """
    taxid: str
    """
    Taxonomy identifier
    """
    uniprot_mappings: list[UniprotMapping]
