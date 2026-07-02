from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal

from cyclopts import Parameter
from cyclopts.types import PositiveInt
from cyclopts.validators import Number
from distributed.deploy.cluster import Cluster

from protein_quest.utils import CopyMethod


@Parameter(name="*")
@dataclass(frozen=True, slots=True)
class CombinedFilterQuery:
    """Query object to apply combined filtering.

    Parameters:
        confidence: Minimal confidence (plDDT) for AlphaFold structures to pass the filter.
        min_residues: Min residues in chain A.
        max_residues: Max residues in chain A.
        minimal_geometry_quality: Minimal geometry quality score to pass the filter.
        top_uniprot_cluster: Maximum number of files to keep for structures per cluster per Uniprot accession.
            Alphafold structures are excluded from this limit.
            The top N structures in each cluster of the resolution clusters.
            The top N structures in each cluster of the PDBe quality clusters.
        top_non_uniprot: Maximum number of files to keep for structures without Uniprot accession.

    """

    confidence: Annotated[float, Parameter(validator=(Number(lte=100, gte=0)))] = (70.0,)
    min_residues: PositiveInt = (0,)
    max_residues: PositiveInt = (10_000_000,)
    minimal_geometry_quality: float = (50.0,)
    top_uniprot_cluster: PositiveInt | None = (1_000,)
    top_non_uniprot: PositiveInt | None = (0,)


def combined_filter(
    input_files: list[Path],
    quality_json: Path,  # TODO change to list of dataclass objects
    output_dir: Path,
    filters: CombinedFilterQuery,
    copy_method: CopyMethod = "copy",
    scheduler_address: str | Cluster | Literal["sequential"] | None = None,
):
    """

    Rules:
    * All non-AlphaFold structures are filtered by number of residues in chain A.
        See `protein-quest filter residue --help` for details.
    * AlphaFold structures are filtered by confidence (plDDT) and afterwards by number of residues in chain A.
        See `protein-quest filter confidence --help` for details.
    * Structures with uniprot accession and resolution are filtered
        by grouping/clustering and  sorting cluster members by resolution.
        See `protein-quest filter resolution --help` for details.
    * Structures with Uniprot accesion and without resolution are filtered by grouping/clustering
        and sorting cluster members by PDBe quality scores.
    * Structures without Uniprot accession and without resolution are filtered by PDBe quality scores.
    * Structures without Uniprot accession and with resolution are filtered/sorted by resolution.

    ```mermaid
    flowchart TD
        A[Input PDB/mmCIF files] --> B{AlphaFold structure?}

        B -->|Yes| C[Filter by confidence plDDT]
        C --> D[Filter by residues in chain A]

        B -->|No| F[Filter by residues in chain A]
        F --> G{UniProt accession?}

        G -->|Yes| H{Resolution available?}
        H -->|Yes| I[Group by UniProt accession and cluster by residue ranges]
        I --> J[Sort cluster members by resolution]
        J --> K[Keep up to top_uniprot_cluster per cluster]

        H -->|No| L[Group by UniProt accession and cluster by residue ranges]
        L --> M[Sort cluster members by PDBe quality]
        M --> N[Keep up to top_uniprot_cluster per cluster]

        G -->|No| O{Resolution available?}
        O -->|No| P[Sort by PDBe quality]
        P --> R[Select up to top_non_uniprot entries]
        O -->|Yes| S[Sort by resolution]
        S --> T[Select up to top_non_uniprot entries]

        D --> Q[Write output structure files]
        K --> Q
        N --> Q
        R --> Q
        T --> Q
    ```

    """
    pass
