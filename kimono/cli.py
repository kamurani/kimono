"""Command line interface for kimono."""

from kimono import __version__
import click as ck 

import pathlib


@ck.command()
@ck.version_option(__version__)
@ck.option(
    "-c",
    "--config_path",
    help="The name the of .yml input config",
    type=ck.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)
def main(config_path):
    """Run kimono."""
    config = parse_config(config_path) if config_path else None 

    pass 
