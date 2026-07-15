import logging

import pytest

from protein_quest.pdbe.clustering import (
    filter_pdbs_on_clustered_resolution,
)
from protein_quest.pdbe.result import PdbResult


def make_pdb(
    pdb_id: str,
    uniprot_chains: str,
    resolution: float | None = 5.0,
    method: str = "X-Ray_Crystallography",
    uniprot_accession: str = "P00000",
) -> PdbResult:
    return PdbResult(
        id=pdb_id,
        uniprot_accession=uniprot_accession,
        method=method,
        resolution=str(resolution) if resolution is not None else None,
        uniprot_chains=uniprot_chains,
    )


@pytest.mark.parametrize(
    "pdbs, expected_ids, expected_log",
    [
        # First example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", 5.4),
                make_pdb("3CCC", "A=1-250", 2.1),
                make_pdb("4DDD", "A=300-400", 8.1),
                make_pdb("5EEE", "A=300-400", 4.6),
                make_pdb("6FFF", "A=500-1000", 1.3),
                make_pdb("7GGG", "A=500-1000", 1.4),
                make_pdb("8HHH", "A=500-1000", 1.6),
            ],
            ["6FFF", "3CCC", "5EEE", "7GGG", "1AAA", "4DDD", "8HHH", "2BBB"],
            [],
            id="three_domains_split",
        ),
        # Second example from https://github.com/haddocking/protein-quest/issues/102
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("4DDD", "A=200-400", 8.1),
                make_pdb("6FFF", "A=500-1000", 1.3),
                make_pdb("9III", "A=1-600", 4.2),
                make_pdb("10JJJ", "A=1-1000", 1.4),
            ],
            ["10JJJ", "6FFF", "1AAA", "9III", "4DDD"],
            [],
            id="overlap_merges",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=-", 1.4),
            ],
            ["1AAA", "2BBB"],
            ["PDB 2BBB is invalid, placing last: Could not determine chain length"],
            id="invalid_chains_last",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=1-250", 3.6),
                make_pdb("2BBB", "A=1-250", resolution=None, method="NMR_Spectroscopy"),
            ],
            ["1AAA", "2BBB"],
            ["PDB 2BBB is invalid, placing last: Resolution is unset"],
            id="unset_resolution_last",
        ),
        pytest.param(
            [
                make_pdb("1AAA", "A=-", 1.4),
                make_pdb("2BBB", "A=1-250", resolution=None, method="NMR_Spectroscopy"),
            ],
            ["2BBB", "1AAA"],
            [
                "PDB 2BBB is invalid, placing last: Resolution is unset",
                "PDB 1AAA is invalid, placing last: Could not determine chain length",
            ],
            id="unset_resolution_last",
        ),
    ],
)
def test_filter_pdbs_on_clustered_resolution(
    pdbs: list[PdbResult], expected_ids: list[str], expected_log: list[str], caplog: pytest.LogCaptureFixture
):
    caplog.set_level(logging.INFO)

    filtered_pdbs = filter_pdbs_on_clustered_resolution(pdbs, top=len(expected_ids))

    filtered_ids = [pdb.id for pdb in filtered_pdbs]
    assert filtered_ids == expected_ids
    for log_entry in expected_log:
        assert log_entry in caplog.text


def test_filter_pdbs_on_clustered_resolution_negative_top():
    with pytest.raises(ValueError, match="Top must be a positive integer"):
        filter_pdbs_on_clustered_resolution([], top=-42)
