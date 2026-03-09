import pytest
from cattrs import ClassValidationError

from protein_quest.converter import converter
from protein_quest.pdbe_3dbeacons.fetch import (
    PruneOptions,
    flatten_structure_summaries,
    search_structure_provider_choices,
    uniprots2structures,
)


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_uniprots2structures_single_provider():
    raw_summaries = await uniprots2structures({"P05067"}, PruneOptions(providers={"alphafill"}))

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P05067",
            "model_url": "https://alphafill.eu/v1/aff/P05067",
            "residue_count": 770,
            "provider": "alphafill",
            "uniprot_accession": "P05067",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.vcr
async def test_uniprots2structures_all_providers():
    raw_summaries = await uniprots2structures(
        {"P38634"},
        PruneOptions(providers=search_structure_provider_choices),
    )

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "C",
            "model_format": "MMCIF",
            "model_identifier": "6g86",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/6g86_updated.cif",
            "residue_count": 16,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "B",
            "model_format": "MMCIF",
            "model_identifier": "3v7d",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/3v7d_updated.cif",
            "residue_count": 19,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "X",
            "model_format": "MMCIF",
            "model_identifier": "8k0g",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/8k0g_updated.cif",
            "residue_count": 4,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "X",
            "model_format": "MMCIF",
            "model_identifier": "8k0g",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/8k0g_updated.cif",
            "residue_count": 43,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e002/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e003",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e003/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00014e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00014/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00014e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00014/ensembles/e002/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00014e003",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00014/ensembles/e003/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00023e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00023/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00023e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00023/ensembles/e002/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00023e003",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00023/ensembles/e003/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00159e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00159/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00160e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00160/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00161e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00161/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00423e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00423/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00424e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00424/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00454e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00454/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00455e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00455/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00486e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00486/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00487e002",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00487/ensembles/e002/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "residue_count": 284,
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "residue_count": 284,
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.default_cassette("test_uniprots2structures_all_providers.yaml")
@pytest.mark.vcr
async def test_uniprots2structures_all_providers_limit1():
    raw_summaries = await uniprots2structures(
        {"P38634"},
        PruneOptions(providers=search_structure_provider_choices, limit=1),
    )

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "C",
            "model_format": "MMCIF",
            "model_identifier": "6g86",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/6g86_updated.cif",
            "residue_count": 16,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "residue_count": 284,
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "residue_count": 284,
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
    raw_summaries = await uniprots2structures({"P38634"}, PruneOptions(providers=providers, limit=1))

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "PED00001e001",
            "model_url": "https://deposition.proteinensemble.org/api/v1/entries/PED00001/ensembles/e001/ensemble-sample/",
            "residue_count": 90,
            "provider": "ped",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "residue_count": 284,
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "residue_count": 284,
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.default_cassette("test_uniprots2structures_all_providers.yaml")
@pytest.mark.vcr
async def test_uniprots2structures_min_residues():
    raw_summaries = await uniprots2structures(
        {"P38634"},
        PruneOptions(providers=search_structure_provider_choices, min_residues=200),
    )

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "residue_count": 284,
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "residue_count": 284,
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.default_cassette("test_uniprots2structures_all_providers.yaml")
@pytest.mark.vcr
async def test_uniprots2structures_max_residues():
    raw_summaries = await uniprots2structures(
        {"P38634"},
        PruneOptions(providers=search_structure_provider_choices, max_residues=16),
    )

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "C",
            "model_format": "MMCIF",
            "model_identifier": "6g86",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/6g86_updated.cif",
            "residue_count": 16,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "X",
            "model_format": "MMCIF",
            "model_identifier": "8k0g",
            "model_url": "https://www.ebi.ac.uk/pdbe/static/entry/8k0g_updated.cif",
            "residue_count": 4,
            "provider": "pdbe",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.default_cassette("test_uniprots2structures_all_providers.yaml")
@pytest.mark.vcr
async def test_uniprots2structures_minandmax_residues():
    raw_summaries = await uniprots2structures(
        {"P38634"},
        PruneOptions(providers=search_structure_provider_choices, min_residues=200, max_residues=300),
    )

    summaries = list(flatten_structure_summaries(raw_summaries))

    expected = [
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "AF-P38634-F1",
            "model_url": "https://alphafold.ebi.ac.uk/files/AF-P38634-F1-model_v6.cif",
            "residue_count": 284,
            "provider": "alphafold",
            "uniprot_accession": "P38634",
        },
        {
            "chain": "A",
            "model_format": "MMCIF",
            "model_identifier": "P38634",
            "model_url": "https://alphafill.eu/v1/aff/P38634",
            "residue_count": 284,
            "provider": "alphafill",
            "uniprot_accession": "P38634",
        },
    ]
    assert summaries == expected


@pytest.mark.asyncio
@pytest.mark.parametrize("batch_size", [0, 1001])
async def test_uniprots2structures_invalid_batch_size(batch_size):
    with pytest.raises(ValueError, match=f"Batch size {batch_size} must be between 1 and 10."):
        await uniprots2structures({"P05067"}, PruneOptions(providers={"alphafill"}), batch_size=batch_size)


class TestPruneOptionsViaCattrsStructure:
    def test_defaults_valid(self):
        options = converter.structure({}, PruneOptions)

        expected = PruneOptions(
            providers={"pdbe", "alphafold"},
            limit=10_000,
            min_residues=None,
            max_residues=None,
        )
        assert options == expected

    def test_all_valid(self):
        options = converter.structure(
            {
                "providers": {"pdbe"},
                "limit": 100,
                "min_residues": 50,
                "max_residues": 500,
            },
            PruneOptions,
        )

        expected = PruneOptions(
            providers={"pdbe"},
            limit=100,
            min_residues=50,
            max_residues=500,
        )
        assert options == expected

    def test_providers_empty(self):
        with pytest.raises(ClassValidationError):
            converter.structure(
                {"providers": []},
                PruneOptions,
            )

    def test_providers_invalid_value(self):
        with pytest.raises(ClassValidationError):
            converter.structure(
                {"providers": ["invalid-provider"]},
                PruneOptions,
            )

    def test_providers_invalid_type(self):
        with pytest.raises(ClassValidationError):
            converter.structure(
                {"providers": [42]},
                PruneOptions,
            )

    def test_limit_negative(self):
        with pytest.raises(ClassValidationError):
            converter.structure(
                {"providers": ["pdbe"], "limit": -1},
                PruneOptions,
            )

    def test_min_residues_negative(self):
        with pytest.raises(ClassValidationError):
            converter.structure(
                {"providers": ["pdbe"], "min_residues": -1},
                PruneOptions,
            )

    def test_max_residues_negative(self):
        with pytest.raises(ClassValidationError):
            converter.structure(
                {"providers": ["pdbe"], "max_residues": -1},
                PruneOptions,
            )
