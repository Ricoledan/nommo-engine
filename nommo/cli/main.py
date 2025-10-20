"""
Main CLI entry point for Nommo Engine.

Provides commands for creating, running, and analyzing simulations.
"""

from pathlib import Path

import click
from rich.console import Console

from nommo.cli import analyze, run, universe
from nommo.core.types import SimulationConfig

console = Console()


@click.group()
@click.version_option(version="0.1.0")
@click.option("--config", type=click.Path(exists=True), help="Configuration file path")
@click.option("--data-dir", type=click.Path(), default="data", help="Data directory path")
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    default="INFO",
    help="Logging level",
)
@click.pass_context
def cli(ctx, config, data_dir, log_level):
    """
    Nommo Engine - A scientific life emergence simulator.

    Explore how self-replicating systems emerge from non-living
    chemistry using real physics and thermodynamics.
    """
    from nommo.utils.logging import setup_logging

    setup_logging(level=log_level)

    ctx.ensure_object(dict)
    ctx.obj["config"] = SimulationConfig(data_dir=Path(data_dir), log_level=log_level)

    if config:
        pass


cli.add_command(universe.universe)
cli.add_command(run.run_cmd)
cli.add_command(analyze.analyze)


if __name__ == "__main__":
    cli()
