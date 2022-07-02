"""Command line interface for CASM."""

import pathlib
import click as c

from CASM.utils.config_parser import parse_config

@c.command
@c.option(
    "-c",
    "--config-path",
    help="The path to the `.yml` config file",
    type=c.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)

def main(config_path):
    """Run job"""
    config = parse_config(path=config_path) if config_path else None 
    c.echo("Running...")
    c.echo(config)




if __name__ == "__main__":
    main()

