"""Command line interface for CASM."""

import pathlib
import click

from CASM.utils.config_parser import parse_config

@click.command
@click.option(
    "-c",
    "--config-path",
    help="The path to the `.yml` config file",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, path_type=pathlib.Path
    ),
)

def main(config_path):
    """Run job"""
    config = parse_config(path=config_path) if config_path else None 
    click.echo("Running...")
    print(config)



if __name__ == "__main__":
    main()

