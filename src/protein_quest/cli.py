import argparse
import logging

from rich.logging import RichHandler

from protein_quest.__version__ import __version__


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Protein Quest CLI", prog="protein-quest")
    parser.add_argument("--log-level", default="WARNING", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)
    return parser



def main():
    parser = make_parser()

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level, handlers=[RichHandler(show_level=False)])
