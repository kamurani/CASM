"""Command line interface for CASM."""

import click

@click.command
def main():
    """Run app"""
    click.echo("CASM")



if __name__ == "__main__":
    main()

