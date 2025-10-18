"""
Simulation execution commands.
"""

import click
import time
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
from rich.live import Live
from rich.panel import Panel
from rich.layout import Layout
from rich.table import Table

from nommo.core.universe import Universe
from nommo.storage.persistence import HDF5Storage, CheckpointManager

console = Console()


@click.group(name='run')
def run_cmd():
    """Execute simulations."""
    pass


@run_cmd.command()
@click.argument('universe_id')
@click.option('--ticks', type=int, default=1000, help='Number of steps')
@click.option('--watch', is_flag=True, help='Live progress display')
@click.option('--checkpoint-interval', type=int, default=1000, help='Steps between checkpoints')
@click.pass_context
def start(ctx, universe_id, ticks, watch, checkpoint_interval):
    """
    Run a simulation.
    
    Examples:
        nommo run start earth-1 --ticks 10000 --watch
        nommo run start my-universe --ticks 5000
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
    universe = storage.load_universe(filepath, found['universe_id'])
    
    checkpoint_mgr = CheckpointManager(storage, checkpoint_interval)
    
    console.print(f"[bold cyan]Starting simulation[/bold cyan]")
    console.print(f"Universe: {universe.params.name}")
    console.print(f"Ticks: {ticks}")
    console.print()
    
    if watch:
        _run_with_display(universe, ticks, checkpoint_mgr)
    else:
        _run_simple(universe, ticks, checkpoint_mgr)
    
    filepath = storage.save_universe(universe)
    console.print(f"\n✅ Simulation complete. Saved to {filepath}")


def _run_simple(universe: Universe, ticks: int, checkpoint_mgr: CheckpointManager):
    """Run simulation with simple progress bar."""
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
    ) as progress:
        task = progress.add_task(f"Running simulation", total=ticks)
        
        start_time = time.perf_counter()
        
        for i in range(ticks):
            metrics = universe.step()
            
            if i % 10 == 0:
                progress.update(task, advance=10)
            
            if i % 100 == 0:
                checkpoint_mgr.checkpoint(universe)
                elapsed = time.perf_counter() - start_time
                tps = i / elapsed if elapsed > 0 else 0
                progress.update(
                    task,
                    description=f"Running ({tps:.1f} ticks/s)"
                )


def _run_with_display(universe: Universe, ticks: int, checkpoint_mgr: CheckpointManager):
    """Run simulation with live display."""
    layout = Layout()
    layout.split_column(
        Layout(name="header", size=3),
        Layout(name="main"),
        Layout(name="footer", size=3)
    )
    
    with Live(layout, refresh_per_second=4) as live:
        start_time = time.perf_counter()
        
        for i in range(ticks):
            metrics = universe.step()
            
            if i % 10 == 0:
                elapsed = time.perf_counter() - start_time
                tps = i / elapsed if elapsed > 0 else 0
                
                layout["header"].update(
                    Panel(
                        f"[bold cyan]Nommo Engine[/bold cyan] - {universe.params.name}",
                        style="cyan"
                    )
                )
                
                stats_table = Table(show_header=False, box=None, padding=(0, 2))
                stats_table.add_column()
                stats_table.add_column()
                stats_table.add_column()
                stats_table.add_column()
                
                stats_table.add_row(
                    f"[cyan]Tick:[/cyan] {universe.tick}",
                    f"[cyan]Time:[/cyan] {universe.time:.2f} ps",
                    f"[cyan]Speed:[/cyan] {tps:.1f} ticks/s",
                    f"[cyan]Progress:[/cyan] {i}/{ticks}"
                )
                
                stats_table.add_row(
                    f"[green]Particles:[/green] {metrics.total_particles}",
                    f"[green]Temperature:[/green] {metrics.current_temperature:.1f} K",
                    f"[green]Energy:[/green] {metrics.total_energy:.1f} kJ/mol",
                    f"[green]Entropy:[/green] {metrics.shannon_entropy:.3f}"
                )
                
                layout["main"].update(Panel(stats_table, title="Metrics"))
                
                layout["footer"].update(
                    Panel(f"Press Ctrl+C to stop", style="dim")
                )
            
            if i % 100 == 0:
                checkpoint_mgr.checkpoint(universe)