"""MkDocs hooks for generating CLI documentation."""

import argparse
import io
import logging

from mkdocs.config.defaults import MkDocsConfig
from mkdocs.structure.files import File, Files

from protein_quest.cli import make_parser

logger = logging.getLogger("mkdocs.plugins.argparse")


def capture_help(parser: argparse.ArgumentParser) -> str:
    """Capture the help text of an argparse parser."""
    with io.StringIO() as help_output:
        parser.print_help(help_output)
        return help_output.getvalue()


def argparser_to_markdown(parser: argparse.ArgumentParser, heading="CLI Reference") -> str:
    prog = parser.prog

    main_help = capture_help(parser)

    lines = [
        f"# {heading}",
        f"Documentation for the `{prog}` script.",
        "```console",
        f"$ {prog} --help",
        main_help,
        "```",
    ]

    subparsers_actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    current_subparsers_action = subparsers_actions[0]

    for sub_cmd_name, sub_cmd_parser in current_subparsers_action.choices.items():
        sub_cmd_help_text = capture_help(sub_cmd_parser)

        lines.extend(
            [
                f"## {sub_cmd_name}",
                "```console",
                f"$ {prog} {sub_cmd_name} --help",
                sub_cmd_help_text,
                "```",
            ]
        )

        # Check for sub-sub-commands
        sub_subparsers_actions = [
            action for action in sub_cmd_parser._actions if isinstance(action, argparse._SubParsersAction)
        ]
        if sub_subparsers_actions:
            sub_current_subparsers_action = sub_subparsers_actions[0]
            for sub_sub_cmd_name, sub_sub_cmd_parser in sub_current_subparsers_action.choices.items():
                sub_sub_cmd_help_text = capture_help(sub_sub_cmd_parser)

                lines.extend(
                    [
                        f"## {sub_cmd_name} {sub_sub_cmd_name}",
                        "```console",
                        f"$ {prog} {sub_cmd_name} {sub_sub_cmd_name} --help",
                        sub_sub_cmd_help_text,
                        "```",
                    ]
                )

    return "\n".join(lines)


def generate_cli_docs() -> str:
    """Generate CLI documentation markdown."""
    parser = make_parser()
    return argparser_to_markdown(parser)


def on_files(files: Files, config: MkDocsConfig) -> Files:
    logger.info("Generating CLI documentation...")
    docs_content = generate_cli_docs()
    cli_md_file = File.generated(
        config=config,
        src_uri="cli.md",
        content=docs_content,
    )
    files.append(cli_md_file)
    return files
