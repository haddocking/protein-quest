import pytest

from protein_quest.taxonomy import Taxon, search_taxon


@pytest.fixture
def expected_human() -> Taxon:
    return Taxon(
        taxon_id="9606",
        scientific_name="Homo sapiens",
        common_name="Human",
        rank="species",
        other_names={
            "Homo sapiens Linnaeus, 1758",
            "human",
            "Home sapiens",
            "Homo sampiens",
            "Homo sapeins",
            "Homo sapian",
            "Homo sapians",
            "Homo sapien",
            "Homo sapience",
            "Homo sapiense",
            "Homo sapients",
            "Homo sapines",
            "Homo spaiens",
            "Homo spiens",
            "Humo sapiens",
            "Homo sapiens (SIRT6)",
            "Homo sapiens (PARIS)",
        },
    )


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_search_taxon(expected_human: Taxon):
    results = await search_taxon("Human", limit=250)

    assert len(results) == 250
    expected0 = expected_human

    results0 = results[0]
    assert results0 == expected0


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_search_taxon_by_id(expected_human: Taxon):
    results = await search_taxon("9606", field="tax_id")

    assert results == [expected_human]
