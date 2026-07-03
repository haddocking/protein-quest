"""Cluster 3D-beacons UniProt summaries by PDBe residue-range overlap.

PDBe Overviews are clustered and pruned with the source-agnostic
[filter_structures_on_clustered_resolution][protein_quest.clustering.filter_structures_on_clustered_resolution].
Overviews from other providers (AlphaFold, SWISS-MODEL, ...) pass through
unchanged.

Attributes:
    PDBE_PROVIDER_RESPONSE: Provider name used in the Overview summaries to identify entries provided by [PDBe](https://www.ebi.ac.uk/pdbe/).
"""

from dataclasses import dataclass

from protein_quest.clustering import filter_structures_on_clustered_resolution
from protein_quest.pdbe_3dbeacons.model import Overview, UniprotSummary

PDBE_PROVIDER_RESPONSE = "PDBe"


@dataclass(frozen=True, eq=False)
class OverviewClusterableEntry:
    """Adapter exposing the [ClusterableStructure][protein_quest.clustering.ClusterableStructure] protocol.

    Wraps an [Overview][protein_quest.pdbe_3dbeacons.model.Overview] so the
    generic clustering algorithm can consume it.
    """

    id: str
    uniprot_start: int
    uniprot_end: int
    resolution_value: float
    sequence_identity: float
    chain_length: int
    overview: Overview

    @classmethod
    def from_overview(cls, overview: Overview) -> "OverviewClusterableEntry":
        summary = overview.summary
        return cls(
            id=summary.model_identifier,
            uniprot_start=summary.uniprot_start,
            uniprot_end=summary.uniprot_end,
            resolution_value=summary.resolution if summary.resolution is not None else 0.0,
            sequence_identity=summary.sequence_identity,
            chain_length=summary.uniprot_end - summary.uniprot_start + 1,
            overview=overview,
        )


def _filter_pdbe_overviews(overviews: list[Overview], top: int) -> list[Overview]:
    entries = [OverviewClusterableEntry.from_overview(o) for o in overviews]
    selected = filter_structures_on_clustered_resolution(entries, top=top)
    return [entry.overview for entry in selected]


def cluster_overviews_per_uniprot(summaries: list[UniprotSummary], top: int) -> list[UniprotSummary]:
    """Cluster PDBe Overviews of each UniProt summary and keep at most ``top`` per accession cluster.

    Overviews from providers other than PDBe are passed through unchanged
    and appended after the pruned PDBe Overviews.

    Args:
        summaries: UniProt summaries to prune.
        top: Maximum number of PDBe Overviews to retain per UniProt accession cluster.

    Returns:
        New list of [UniprotSummary][protein_quest.pdbe_3dbeacons.model.UniprotSummary] with PDBe Overviews pruned.

    Raises:
        ValueError: If ``top`` is not a positive integer.
    """
    if top <= 0:
        msg = "Top must be a positive integer."
        raise ValueError(msg)

    pruned: list[UniprotSummary] = []
    for summary in summaries:
        if summary.structures is None:
            pruned.append(summary)
            continue
        pdbe_overviews = [o for o in summary.structures if o.summary.provider == PDBE_PROVIDER_RESPONSE]
        other_overviews = [o for o in summary.structures if o.summary.provider != PDBE_PROVIDER_RESPONSE]
        if pdbe_overviews:
            pdbe_overviews = _filter_pdbe_overviews(pdbe_overviews, top=top)
        pruned.append(
            UniprotSummary(
                uniprot_entry=summary.uniprot_entry,
                structures=pdbe_overviews + other_overviews,
            )
        )
    return pruned
