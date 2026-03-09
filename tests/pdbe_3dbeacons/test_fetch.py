from typing import cast

import pytest

from protein_quest.pdbe_3dbeacons.fetch import (
    Provider,
    flatten_structure_summaries,
    search_structure_provider_choices,
    uniprots2structures,
)


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_uniprots2structures_single_provider():
    raw_summaries = await uniprots2structures({"P05067"}, {"alphafill"})

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_uniprots2structures_all_providers():
    raw_summaries = await uniprots2structures({"P38634"}, search_structure_provider_choices)

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "C",
            "model_format": "MMCIF",
            "model_identifier": "6g86",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/6g86_updated.cif",
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "B",
            "model_format": "MMCIF",
            "model_identifier": "3v7d",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/3v7d_updated.cif",
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "X",
            "model_format": "MMCIF",
            "model_identifier": "8k0g",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/8k0g_updated.cif",
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "X",
            "model_format": "MMCIF",
            "model_identifier": "8k0g",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/8k0g_updated.cif",
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e002/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e003",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e003/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00014e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00014/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00014e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00014/ensembles/e002/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00014e003",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00014/ensembles/e003/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00023e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00023/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00023e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00023/ensembles/e002/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00023e003",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00023/ensembles/e003/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00159e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00159/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00160e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00160/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00161e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00161/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00423e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00423/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00424e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00424/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00454e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00454/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00455e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00455/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00486e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00486/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00487e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00487/ensembles/e002/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.default_cassette("test_uniprots2structures_all_providers.yaml")
@pytest.mark.vcr
async def test_uniprots2structures_all_providers_limit1():
    raw_summaries = await uniprots2structures({"P38634"}, search_structure_provider_choices, limit=1)

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "C",
            "model_format": "MMCIF",
            "model_identifier": "6g86",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/6g86_updated.cif",
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_uniprots2structures_1providermissing():
    providers = search_structure_provider_choices.copy()
    providers.remove("pdbe")
    raw_summaries = await uniprots2structures({"P38634"}, providers, limit=1)

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e001/ensemble-sample/",
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.parametrize("batch_size", [0, 1001])
async def test_uniprots2structures_invalid_batch_size(batch_size):
    with pytest.raises(ValueError, match=f"Batch size {batch_size} must be between 1 and 10."):
        await uniprots2structures({"P05067"}, {"alphafill"}, batch_size=batch_size)


@pytest.mark.asyncio
async def test_uniprots2structures_invalid_providers():
    providers = cast("set[Provider]", {"invalid_provider"})
    with pytest.raises(
        ValueError, match=f"Providers {providers} must be subset of {search_structure_provider_choices}"
    ):
        await uniprots2structures({"P05067"}, providers)
