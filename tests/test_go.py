import pytest

from protein_quest.go import GoTerm, SearchResponse, converter, search_gene_ontology_term


def test_converter():
    input_data = b'{"numberOfHits":4,"results":[{"id":"GO:0006816","isObsolete":false,"name":"calcium ion transport","definition":{"text":"The directed movement of calcium (Ca) ions into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore."},"aspect":"biological_process"},{"id":"GO:0000811","isObsolete":false,"name":"GINS complex","definition":{"text":"A heterotetrameric protein complex that associates with replication origins, where it is required for the initiation of DNA replication, and with replication forks."},"aspect":"cellular_component"},{"id":"GO:0070966","isObsolete":false,"name":"nuclear-transcribed mRNA catabolic process, no-go decay","definition":{"text":"The chemical reactions and pathways resulting in the breakdown of the transcript body of a nuclear-transcribed mRNA with stalls in translation elongation."},"aspect":"biological_process"},{"id":"GO:0001591","isObsolete":false,"name":"dopamine neurotransmitter receptor activity, coupled via Gi/Go","definition":{"text":"Combining with the neurotransmitter dopamine and activating adenylate cyclase via coupling to Gi/Go to initiate a change in cell activity."},"aspect":"molecular_function"}],"pageInfo":{"resultsPerPage":25,"current":1,"total":1}}'

    result = converter.loads(input_data, SearchResponse)

    assert isinstance(result, SearchResponse)
    assert len(result.results) == 4
    assert result.number_of_hits == 4
    assert (
        result.results[0].definition
        == "The directed movement of calcium (Ca) ions into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore."
    )
    assert result.results[0].is_obsolete is False


@pytest.mark.vcr
@pytest.mark.asyncio
async def test_search_gene_ontology_term():
    results = await search_gene_ontology_term("GO:0043293", limit=1)

    expected = [
        GoTerm(
            id="GO:0043293",
            is_obsolete=False,
            name="apoptosome",
            definition="A multisubunit protein complex involved in the signaling "
            "phase of the apoptotic process. In mammals it is typically "
            "composed of seven Apaf-1 subunits bound to cytochrome c "
            "and caspase-9. A similar complex to promote apoptosis is "
            "formed from homologous gene products in other eukaryotic "
            "organisms.",
            aspect="cellular_component",
        )
    ]
    assert results == expected
