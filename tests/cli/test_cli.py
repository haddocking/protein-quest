import pytest

from protein_quest.cli import main


def test_app_help(capsys: pytest.CaptureFixture[str]):
    # Smoke test to ensure the CLI can be invoked and shows help message
    main(["--help"])

    captured = capsys.readouterr()
    assert "Protein Quest CLI" in captured.out
