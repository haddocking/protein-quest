"""MkDocs hooks for generating CLI documentation.

To get colored output the help should have rich text using rich or rich-argparse package.
"""
# TODO convert into hook plugin in its own repository
# see https://www.mkdocs.org/dev-guide/plugins/#developing-plugins

import argparse
import io
import logging
from textwrap import dedent

from mkdocs.config.defaults import MkDocsConfig
from mkdocs.structure.files import File, Files
from rich.console import Console
from rich.text import Text
from rich_argparse import ArgumentDefaultsRichHelpFormatter

from protein_quest.cli import make_parser

logger = logging.getLogger("mkdocs.plugins.argparse")


def capture_help(parser: argparse.ArgumentParser) -> str:
    """Capture the help text of an argparse parser as HTML."""
    # Based on https://github.com/hamdanal/rich-argparse/blob/e28584ac56ddd46f4079d037c27f24f0ec4eccb4/rich_argparse/_argparse.py#L545
    # but with export instead of save

    # Overwrite default colors as on mkdocs black text is not visible in dark mode
    ArgumentDefaultsRichHelpFormatter.styles["argparse.help"] = "green"
    ArgumentDefaultsRichHelpFormatter.styles["argparse.text"] = "green"

    text = Text.from_ansi(parser.format_help())
    console = Console(file=io.StringIO(), record=True)
    console.print(text, crop=False)
    code_format = dedent("""\
        <pre style="font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">
        <code style="font-family:inherit">
        {code}
        </code>
        </pre>
    """)
    return console.export_html(code_format=code_format, inline_styles=True)


def argparser_to_markdown(parser: argparse.ArgumentParser, heading="CLI Reference") -> str:
    prog = parser.prog

    main_help = capture_help(parser)

    lines = [
        f"# {heading}",
        f"Documentation for the `{prog}` script.",
        "```console",
        f"{prog} --help",
        "```",
        main_help,
    ]

    subparsers_actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
    current_subparsers_action = subparsers_actions[0]

    for sub_cmd_name, sub_cmd_parser in current_subparsers_action.choices.items():
        sub_cmd_help_text = capture_help(sub_cmd_parser)

        lines.extend(
            [
                f"## {sub_cmd_name}",
                "```console",
                f"{prog} {sub_cmd_name} --help",
                "```",
                sub_cmd_help_text,
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
                        f"{prog} {sub_cmd_name} {sub_sub_cmd_name} --help",
                        "```",
                        sub_sub_cmd_help_text,
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
