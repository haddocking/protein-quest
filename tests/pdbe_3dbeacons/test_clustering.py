"""Tests for clustering 3D-beacons UniProt summaries."""

import pytest

from protein_quest.pdbe_3dbeacons.clustering import OverviewClusterableEntry, cluster_overviews_per_uniprot
from protein_quest.pdbe_3dbeacons.model import (
    AppUniprotSchemaEntity,
    AppUniprotSchemaSummaryItems,
    Overview,
    UniprotEntry,
    UniprotSummary,
)


def _entity() -> AppUniprotSchemaEntity:
    return AppUniprotSchemaEntity(
        entity_type="POLYMER",
        description="protein",
        chain_ids=["A"],
    )


def _pdbe_overview(
    model_identifier: str,
    uniprot_start: int,
    uniprot_end: int,
    resolution: float | None = None,
    sequence_identity: float = 1.0,
) -> Overview:
    return Overview(
        summary=AppUniprotSchemaSummaryItems(
            model_identifier=model_identifier,
            model_category="EXPERIMENTALLY DETERMINED",
            model_url=f"https://example.org/{model_identifier}.cif",
            model_format="MMCIF",
            provider="PDBe",
            created="2024-01-01",
            sequence_identity=sequence_identity,
            uniprot_start=uniprot_start,
            uniprot_end=uniprot_end,
            coverage=1.0,
            entities=[_entity()],
            resolution=resolution,
        )
    )


def _alphafold_overview(model_identifier: str, uniprot_start: int, uniprot_end: int) -> Overview:
    return Overview(
        summary=AppUniprotSchemaSummaryItems(
            model_identifier=model_identifier,
            model_category="TEMPLATE-BASED",
            model_url=f"https://alphafold.ebi.ac.uk/{model_identifier}.cif",
            model_format="MMCIF",
            provider="AlphaFold DB",
            created="2024-01-01",
            sequence_identity=1.0,
            uniprot_start=uniprot_start,
            uniprot_end=uniprot_end,
            coverage=1.0,
            entities=[_entity()],
        )
    )


def test_OverviewClusterableEntry_is_hashable():
    overview = _pdbe_overview("1aaa", 1, 250, resolution=3.6)
    entry = OverviewClusterableEntry.from_overview(overview)
    assert isinstance(hash(entry), int)


def _summary(accession: str, overviews: list[Overview]) -> UniprotSummary:
    return UniprotSummary(uniprot_entry=UniprotEntry(ac=accession), structures=overviews)


def test_cluster_overviews_per_uniprot_keeps_top_pdbe_and_passes_alphafold_through():
    # Two PDBe overviews covering the same UniProt residue range (1-250) -> one cluster.
    # 3CCC has the best (lowest) resolution so should be retained with top=1.
    # 4DDD covers a non-overlapping range -> separate cluster, not selected at top=1.
    # The AlphaFold entry must always be kept regardless of clustering.
    summary = _summary(
        "P12345",
        [
            _pdbe_overview("1aaa", 1, 250, resolution=3.6),
            _pdbe_overview("3ccc", 1, 250, resolution=2.1),
            _pdbe_overview("4ddd", 300, 400, resolution=8.1),
            _alphafold_overview("AF-P12345-F1", 1, 400),
        ],
    )

    [pruned] = cluster_overviews_per_uniprot([summary], top=1)

    assert pruned.uniprot_entry is not None
    assert pruned.uniprot_entry.ac == "P12345"
    assert pruned.structures is not None
    pdbe_ids = [o.summary.model_identifier for o in pruned.structures if o.summary.provider == "PDBe"]
    alphafold_ids = [o.summary.model_identifier for o in pruned.structures if o.summary.provider == "AlphaFold DB"]
    assert pdbe_ids == ["3ccc"]
    assert alphafold_ids == ["AF-P12345-F1"]


def test_cluster_overviews_per_uniprot_passes_through_when_no_pdbe():
    summary = _summary(
        "P99999",
        [
            _alphafold_overview("AF-P99999-F1", 1, 400),
        ],
    )

    [pruned] = cluster_overviews_per_uniprot([summary], top=1)

    assert pruned.structures is not None
    assert [o.summary.model_identifier for o in pruned.structures] == ["AF-P99999-F1"]


def test_cluster_overviews_per_uniprot_handles_none_structures():
    summary = UniprotSummary(uniprot_entry=UniprotEntry(ac="P00000"), structures=None)

    [pruned] = cluster_overviews_per_uniprot([summary], top=1)

    assert pruned.structures is None


def test_cluster_overviews_per_uniprot_rejects_non_positive_top():
    with pytest.raises(ValueError, match="positive integer"):
        cluster_overviews_per_uniprot([], top=0)
