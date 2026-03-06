import pytest

from protein_quest.pdbe_3dbeacons.fetch import flatten, uniprots2structures


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_uniprots2structures():
    raw_summaries = await uniprots2structures({"P05067"}, {"alphafill"})

    summaries = list(flatten(raw_summaries))

    expected = [
        {
            "chains": "A",
            "description": "Amyloid-beta precursor protein",
            "entity_type": "POLYMER",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
        {
            "chains": "B:C:D:E:F:M",
            "description": "ZINC ION",
            "entity_type": "NON-POLYMER",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
        {
            "chains": "G",
            "description": "MAGNESIUM ION",
            "entity_type": "NON-POLYMER",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
        {
            "chains": "H:I:L",
            "description": "COPPER (II) ION",
            "entity_type": "NON-POLYMER",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
        {
            "chains": "J:K",
            "description": "(3R)-butane-1,3-diol",
            "entity_type": "NON-POLYMER",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
        {
            "chains": "N",
            "description": "PROTOPORPHYRIN IX CONTAINING FE",
            "entity_type": "NON-POLYMER",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
    ]
    assert summaries == expected
