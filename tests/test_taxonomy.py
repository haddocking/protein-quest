import pytest

from protein_quest.taxonomy import Taxon, search_taxon


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_search_taxon():
    results = await search_taxon("Human", limit=250)

    assert len(results) == 250
    expected0 = Taxon(
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

    results0 = results[0]
    assert results0 == expected0
