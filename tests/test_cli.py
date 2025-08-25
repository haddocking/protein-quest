import pytest

from protein_quest.cli import make_parser


def test_make_parser_help(capsys: pytest.CaptureFixture[str]):
    in_args = ["--help"]
    parser = make_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(in_args)

    captured = capsys.readouterr()
    assert "Protein Quest CLI" in captured.out
    # TODO like ../docs/cli_doc_hook.py loop over all sub commands to get help
