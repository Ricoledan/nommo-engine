"""
Universe management commands.
"""

import click
import json
import yaml
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.prompt import Prompt, Confirm, IntPrompt, FloatPrompt

from nommo.core.types import UniverseParameters, ParticleTypeConfig
from nommo.core.universe import Universe
from nommo.storage.persistence import HDF5Storage
from nommo.presets import get_preset

console = Console()


@click.group()
def universe():
    """Create and manage simulation universes."""
    pass


@universe.command()
@click.argument('name', required=False)
@click.option('--preset', help='Use a preset configuration')
@click.option('--interactive', '-i', is_flag=True, help='Interactive setup')
@click.option('--config', type=click.Path(exists=True), help='Load from config file')
@click.pass_context
def create(ctx, name, preset, interactive, config):
    """
    Create a new universe.
    
    Examples:
        nommo universe create earth-1 --preset earth-like
        nommo universe create -i
        nommo universe create my-universe --config config.yaml
    """
    storage_dir = ctx.obj['config'].data_dir / 'universes'
    storage = HDF5Storage(storage_dir)
    
    if config:
        with open(config) as f:
            if config.endswith('.yaml') or config.endswith('.yml'):
                data = yaml.safe_load(f)
            else:
                data = json.load(f)
        params = UniverseParameters(**data)
        
    elif interactive:
        console.print("[bold cyan]Creating new universe[/bold cyan]")
        
        name = Prompt.ask("Universe name", default="my-universe")
        description = Prompt.ask("Description", default="")
        
        preset_choice = Prompt.ask(
            "Start from preset?",
            choices=["earth-like", "high-energy", "minimal", "custom"],
            default="earth-like"
        )
        
        if preset_choice != "custom":
            params = get_preset(preset_choice)
            params.name = name
            params.description = description
        else:
            params = _interactive_setup(name, description)
            
    elif preset:
        params = get_preset(preset)
        params.name = name or f"{preset}-universe"
        
    else:
        if not name:
            console.print("[red]Universe name required[/red]")
            return
            
        params = get_preset("minimal")
        params.name = name
    
    uni = Universe(params)
    
    filepath = storage.save_universe(uni)
    
    console.print(f"✅ Universe '{params.name}' created")
    console.print(f"   ID: {uni.id}")
    console.print(f"   Particles: {len(uni.particles)}")
    console.print(f"   Saved to: {filepath}")


@universe.command()
@click.option('--verbose', '-v', is_flag=True, help='Show detailed information')
@click.pass_context
def list(ctx, verbose):
    """
    List all universes.
    
    Examples:
        nommo universe list
        nommo universe list -v
    """
    storage_dir = ctx.obj['config'].data_dir / 'universes'
    storage = HDF5Storage(storage_dir)
    
    saves = storage.list_saves()
    
    if not saves:
        console.print("No universes found")
        return
    
    table = Table(title="Active Universes")
    table.add_column("Name", style="cyan")
    table.add_column("ID", style="magenta")
    table.add_column("Tick", justify="right")
    table.add_column("Time (ps)", justify="right")
    table.add_column("Created", style="green")
    
    for save in saves:
        table.add_row(
            save['filename'].replace('.h5', ''),
            save['universe_id'][:8],
            str(save['tick']),
            f"{save['time']:.2f}",
            save['creation_time'][:19] if save['creation_time'] else "Unknown"
        )
    
    console.print(table)


@universe.command()
@click.argument('universe_id')
@click.pass_context
def show(ctx, universe_id):
    """
    Show universe details.
    
    Examples:
        nommo universe show earth-1
        nommo universe show abc123
    """
    storage_dir = ctx.obj['config'].data_dir / 'universes'
    storage = HDF5Storage(storage_dir)
    
    saves = storage.list_saves()
    
    found = None
    for save in saves:
        if universe_id in save['filename'] or universe_id in save['universe_id']:
            found = save
            break
    
    if not found:
        console.print(f"[red]Universe '{universe_id}' not found[/red]")
        return
    
    filepath = storage_dir / found['filename']
    uni = storage.load_universe(filepath, found['universe_id'])
    
    console.print(f"[bold cyan]Universe: {uni.params.name}[/bold cyan]")
    console.print(f"ID: {uni.id}")
    console.print(f"Description: {uni.params.description or 'N/A'}")
    console.print()
    
    console.print("[bold]Status:[/bold]")
    console.print(f"  Tick: {uni.tick}")
    console.print(f"  Time: {uni.time:.3f} ps")
    console.print(f"  Particles: {len(uni.particles)}")
    console.print()
    
    console.print("[bold]Physics:[/bold]")
    console.print(f"  Temperature: {uni.params.thermodynamics.temperature} K")
    console.print(f"  Box: {uni.params.box.dimensions}")
    console.print(f"  Timestep: {uni.params.physics.timestep} ps")
    console.print()
    
    if uni.metrics:
        console.print("[bold]Metrics:[/bold]")
        console.print(f"  Total Energy: {uni.metrics.total_energy:.2f} kJ/mol")
        console.print(f"  Current Temp: {uni.metrics.current_temperature:.1f} K")
        console.print(f"  Entropy: {uni.metrics.shannon_entropy:.3f}")
        console.print(f"  Bonds: {uni.metrics.total_bonds}")


@universe.command()
@click.argument('universe_id')
@click.option('--force', '-f', is_flag=True, help='Skip confirmation')
@click.pass_context
def delete(ctx, universe_id, force):
    """Delete a universe."""
    storage_dir = ctx.obj['config'].data_dir / 'universes'
    
    saves = HDF5Storage(storage_dir).list_saves()
    
    found = None
    for save in saves:
        if universe_id in save['filename'] or universe_id in save['universe_id']:
            found = save
            break
    
    if not found:
        console.print(f"[red]Universe '{universe_id}' not found[/red]")
        return
    
    if not force:
        if not Confirm.ask(f"Delete universe '{found['filename']}'?"):
            return
    
    filepath = storage_dir / found['filename']
    filepath.unlink()
    
    console.print(f"✅ Universe deleted: {found['filename']}")


def _interactive_setup(name: str, description: str) -> UniverseParameters:
    """Interactive universe parameter setup."""
    console.print("\n[bold]Physics Parameters:[/bold]")
    temperature = FloatPrompt.ask("Temperature (K)", default=300.0)
    timestep = FloatPrompt.ask("Timestep (ps)", default=0.001)
    
    console.print("\n[bold]Box Parameters:[/bold]")
    box_size = FloatPrompt.ask("Box size (nm)", default=10.0)
    
    console.print("\n[bold]Particle Types:[/bold]")
    n_types = IntPrompt.ask("Number of particle types", default=2)
    
    particle_types = {}
    initial_composition = {}
    
    for i in range(n_types):
        console.print(f"\n[cyan]Type {i+1}:[/cyan]")
        ptype_name = Prompt.ask("  Name", default=f"type_{i+1}")
        mass = FloatPrompt.ask("  Mass (amu)", default=50.0)
        radius = FloatPrompt.ask("  Radius (nm)", default=0.2)
        count = IntPrompt.ask("  Initial count", default=100)
        
        particle_types[ptype_name] = ParticleTypeConfig(
            name=ptype_name,
            mass=mass,
            radius=radius
        )
        initial_composition[ptype_name] = count
    
    return UniverseParameters(
        name=name,
        description=description,
        particle_types=particle_types,
        initial_composition=initial_composition
    )