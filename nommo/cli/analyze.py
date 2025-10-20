"""
Analysis commands for simulation results.
"""


import click
import matplotlib.pyplot as plt
from rich.console import Console
from rich.table import Table

from nommo.storage.persistence import HDF5Storage

console = Console()


@click.group()
def analyze():
    """Analyze simulation results."""
    pass


@analyze.command()
@click.argument("universe_id")
@click.pass_context
def stats(ctx, universe_id):
    """
    Show universe statistics.

    Examples:
        nommo analyze stats earth-1
    """
    storage_dir = ctx.obj["config"].data_dir / "universes"
    storage = HDF5Storage(storage_dir)

    saves = storage.list_saves()
    found = None
    for save in saves:
        if universe_id in save["filename"] or universe_id in save["universe_id"]:
            found = save
            break

    if not found:
        console.print(f"[red]Universe '{universe_id}' not found[/red]")
        return

    filepath = storage_dir / found["filename"]
    universe = storage.load_universe(filepath, found["universe_id"])

    metrics = universe.metrics

    table = Table(title=f"Statistics for {universe.params.name}")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="magenta", justify="right")

    table.add_row("Total Particles", str(metrics.total_particles))
    table.add_row("Total Bonds", str(metrics.total_bonds))
    table.add_row("Avg Bonds/Particle", f"{metrics.avg_bonds_per_particle:.2f}")
    table.add_row("Shannon Entropy", f"{metrics.shannon_entropy:.3f}")
    table.add_row("Temperature", f"{metrics.current_temperature:.1f} K")
    table.add_row("Kinetic Energy", f"{metrics.total_kinetic_energy:.2f} kJ/mol")
    table.add_row("Potential Energy", f"{metrics.total_potential_energy:.2f} kJ/mol")
    table.add_row("Total Energy", f"{metrics.total_energy:.2f} kJ/mol")
    table.add_row("Max Generation", str(metrics.generation_max))
    table.add_row("Complexity Trend", metrics.complexity_trend)

    console.print(table)

    if metrics.particles_by_type:
        type_table = Table(title="Particle Distribution")
        type_table.add_column("Type", style="cyan")
        type_table.add_column("Count", style="magenta", justify="right")
        type_table.add_column("Percentage", style="green", justify="right")

        total = metrics.total_particles
        for ptype, count in metrics.particles_by_type.items():
            percentage = (count / total * 100) if total > 0 else 0
            type_table.add_row(ptype, str(count), f"{percentage:.1f}%")

        console.print(type_table)


@analyze.command()
@click.argument("universe_id")
@click.argument(
    "metric", type=click.Choice(["temperature", "energy", "particles", "entropy", "bonds"])
)
@click.option("--output", help="Save plot to file")
@click.option("--show/--no-show", default=True, help="Display plot")
@click.pass_context
def plot(ctx, universe_id, metric, output, show):
    """
    Plot metric over time.

    Examples:
        nommo analyze plot earth-1 energy
        nommo analyze plot earth-1 temperature --output temp.png
    """
    storage_dir = ctx.obj["config"].data_dir / "universes"

    trajectory_file = storage_dir / f"trajectory_{universe_id[:8]}.h5"

    if not trajectory_file.exists():
        console.print(f"[red]No trajectory data found for '{universe_id}'[/red]")
        return

    import h5py

    with h5py.File(trajectory_file, "r") as f:
        traj = f["trajectory"]
        times = traj["times"][:]

        if metric == "temperature":
            data = traj["temperatures"][:]
            ylabel = "Temperature (K)"
            title = "Temperature Evolution"
        elif metric == "energy":
            kinetic = traj["kinetic_energy"][:]
            potential = traj["potential_energy"][:]
            total = traj["total_energy"][:]
            ylabel = "Energy (kJ/mol)"
            title = "Energy Evolution"
        elif metric == "particles":
            data = traj["particle_count"][:]
            ylabel = "Particle Count"
            title = "Particle Population"
        elif metric == "entropy":
            data = traj["entropy"][:]
            ylabel = "Shannon Entropy"
            title = "System Entropy"
        elif metric == "bonds":
            console.print("[yellow]Bond data not yet available in trajectory[/yellow]")
            return

    plt.figure(figsize=(10, 6))

    if metric == "energy":
        plt.plot(times, kinetic, label="Kinetic", alpha=0.7)
        plt.plot(times, potential, label="Potential", alpha=0.7)
        plt.plot(times, total, label="Total", linewidth=2)
        plt.legend()
    else:
        plt.plot(times, data)

    plt.xlabel("Time (ps)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)

    if output:
        plt.savefig(output, dpi=150, bbox_inches="tight")
        console.print(f"✅ Plot saved to {output}")

    if show:
        plt.show()


@analyze.command()
@click.argument("universe_ids", nargs=-1, required=True)
@click.pass_context
def compare(ctx, universe_ids):
    """
    Compare multiple universes.

    Examples:
        nommo analyze compare earth-1 earth-2 earth-3
    """
    storage_dir = ctx.obj["config"].data_dir / "universes"
    storage = HDF5Storage(storage_dir)

    universes = []
    for uid in universe_ids:
        saves = storage.list_saves()
        for save in saves:
            if uid in save["filename"] or uid in save["universe_id"]:
                filepath = storage_dir / save["filename"]
                uni = storage.load_universe(filepath, save["universe_id"])
                universes.append(uni)
                break

    if not universes:
        console.print("[red]No universes found[/red]")
        return

    table = Table(title="Universe Comparison")
    table.add_column("Metric", style="cyan")

    for uni in universes:
        table.add_column(uni.params.name[:15], justify="right")

    metrics_to_compare = [
        ("Tick", lambda u: str(u.tick)),
        ("Time (ps)", lambda u: f"{u.time:.2f}"),
        ("Particles", lambda u: str(u.metrics.total_particles)),
        ("Temperature (K)", lambda u: f"{u.metrics.current_temperature:.1f}"),
        ("Total Energy", lambda u: f"{u.metrics.total_energy:.1f}"),
        ("Entropy", lambda u: f"{u.metrics.shannon_entropy:.3f}"),
        ("Bonds", lambda u: str(u.metrics.total_bonds)),
    ]

    for metric_name, metric_func in metrics_to_compare:
        row = [metric_name]
        for uni in universes:
            try:
                row.append(metric_func(uni))
            except Exception:
                row.append("N/A")
        table.add_row(*row)

    console.print(table)
