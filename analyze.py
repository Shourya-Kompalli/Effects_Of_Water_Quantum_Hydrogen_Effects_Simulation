"""Analyze q-TIP4P/F PIMD distance distributions and ring-polymer spread."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    frames: list[np.ndarray] = []
    symbols: list[str] | None = None

    with path.open(encoding="utf-8") as trajectory:
        while True:
            atom_count_line = trajectory.readline()
            if not atom_count_line:
                break
            if not atom_count_line.strip():
                continue

            atom_count = int(atom_count_line.strip())
            comment = trajectory.readline()
            if not comment:
                raise ValueError(f"Missing XYZ comment line in {path}.")

            frame_symbols: list[str] = []
            positions: list[tuple[float, float, float]] = []
            for _ in range(atom_count):
                fields = trajectory.readline().split()
                if len(fields) < 4:
                    raise ValueError(f"Malformed XYZ atom record in {path}.")
                frame_symbols.append(fields[0])
                positions.append(tuple(float(value) for value in fields[1:4]))

            if symbols is None:
                symbols = frame_symbols
            elif frame_symbols != symbols:
                raise ValueError(f"Atom ordering changes between frames in {path}.")
            frames.append(np.asarray(positions, dtype=float))

    if symbols is None or not frames:
        raise ValueError(f"No XYZ frames found in {path}.")
    return symbols, np.stack(frames)


def read_cubic_cell(pdb_path: Path) -> float:
    with pdb_path.open(encoding="ascii") as structure:
        for line in structure:
            if line.startswith("CRYST1"):
                lengths = [float(line[start : start + 9]) for start in (6, 15, 24)]
                if max(lengths) - min(lengths) > 1e-5:
                    raise ValueError("Analysis currently requires a cubic unit cell.")
                return lengths[0]
    raise ValueError(f"No CRYST1 unit cell found in {pdb_path}.")


def minimum_image(displacements: np.ndarray, cell: float) -> np.ndarray:
    return displacements - cell * np.rint(displacements / cell)


def smooth_density(values: np.ndarray, bins: int = 120) -> tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(0.0, max(float(np.max(values)) * 1.05, 0.1), bins + 1)
    counts, _ = np.histogram(values, bins=edges)
    offsets = np.arange(-5, 6, dtype=float)
    kernel = np.exp(-0.5 * (offsets / 1.5) ** 2)
    kernel /= kernel.sum()
    density = np.convolve(counts.astype(float), kernel, mode="same")
    widths = np.diff(edges)
    density /= np.sum(density * widths)
    return 0.5 * (edges[:-1] + edges[1:]), density


def collect_measurements(
    bead_trajectories: list[np.ndarray],
    symbols: list[str],
    cell: float,
    discard: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if len(symbols) % 3 != 0:
        raise ValueError("Expected O-H-H triplets for every water molecule.")
    for index in range(0, len(symbols), 3):
        if symbols[index : index + 3] != ["O", "H", "H"]:
            raise ValueError("Atoms must be ordered as O-H-H for each water molecule.")

    nframes = bead_trajectories[0].shape[0]
    if discard < 0 or discard >= nframes:
        raise ValueError(f"Discard count must be between 0 and {nframes - 1}.")

    positions = np.stack([trajectory[discard:] for trajectory in bead_trajectories])
    molecule_count = len(symbols) // 3
    oxygen_indices = np.arange(molecule_count) * 3
    hydrogen_indices = np.column_stack(
        (oxygen_indices + 1, oxygen_indices + 2)
    ).ravel()
    hydrogen_molecules = np.repeat(np.arange(molecule_count), 2)
    oxygen_mask = np.asarray(symbols) == "O"
    hydrogen_mask = np.asarray(symbols) == "H"
    bond_distances: list[np.ndarray] = []
    neighbor_distances: list[np.ndarray] = []
    spread_by_element: dict[str, list[np.ndarray]] = {"H": [], "O": []}

    for frame_index in range(positions.shape[1]):
        frame = positions[:, frame_index]
        covalent = frame[:, hydrogen_indices, :] - frame[:, oxygen_indices.repeat(2), :]
        covalent = minimum_image(covalent, cell)
        bond_distances.append(np.linalg.norm(covalent, axis=-1).ravel())

        oxygen_positions = frame[:, oxygen_indices, :]
        hydrogen_positions = frame[:, hydrogen_indices, :]
        intermolecular = (
            oxygen_positions[:, :, None, :] - hydrogen_positions[:, None, :, :]
        )
        intermolecular = minimum_image(intermolecular, cell)
        pair_distances = np.linalg.norm(intermolecular, axis=-1)
        same_molecule = (
            np.arange(molecule_count)[:, None] == hydrogen_molecules[None, :]
        )
        pair_distances[:, same_molecule] = np.inf
        neighbor_distances.append(np.min(pair_distances, axis=1).ravel())

        reference = frame[0]
        relative = minimum_image(frame - reference[None, :, :], cell)
        centroid = np.mean(relative, axis=0)
        bead_offsets = minimum_image(relative - centroid[None, :, :], cell)
        radii = np.linalg.norm(bead_offsets, axis=-1)
        spread_by_element["H"].append(radii[:, hydrogen_mask].ravel())
        spread_by_element["O"].append(radii[:, oxygen_mask].ravel())

    return (
        np.concatenate(bond_distances),
        np.concatenate(neighbor_distances),
        np.concatenate(spread_by_element["H"]),
        np.concatenate(spread_by_element["O"]),
    )


def summarize(name: str, values: np.ndarray) -> dict[str, float | str]:
    return {
        "measurement": name,
        "mean_angstrom": float(np.mean(values)),
        "median_angstrom": float(np.median(values)),
        "p95_angstrom": float(np.percentile(values, 95)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", default="simulation", help="i-PI output prefix.")
    parser.add_argument(
        "--structure",
        type=Path,
        default=Path("data/water_32.pdb"),
        help="Periodic PDB structure used to obtain the fixed cell.",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("analysis"), help="Output directory."
    )
    parser.add_argument(
        "--discard",
        type=int,
        default=0,
        help="Number of initial trajectory frames to omit as equilibration.",
    )
    args = parser.parse_args()

    trajectories = sorted(
        Path(".").glob(f"{args.prefix}.pos_*.xyz"),
        key=lambda path: int(path.stem.rsplit("_", 1)[1]),
    )
    if not trajectories:
        raise FileNotFoundError(
            f"No bead trajectories matching {args.prefix}.pos_*.xyz were found."
        )

    symbols: list[str] | None = None
    bead_positions: list[np.ndarray] = []
    for path in trajectories:
        bead_symbols, positions = read_xyz(path)
        if symbols is None:
            symbols = bead_symbols
        elif bead_symbols != symbols:
            raise ValueError(f"Atom ordering differs in {path}.")
        bead_positions.append(positions)

    if any(positions.shape != bead_positions[0].shape for positions in bead_positions):
        raise ValueError("Bead trajectories must have equal frame and atom counts.")

    cell = read_cubic_cell(args.structure)
    bond_distances, neighbor_distances, hydrogen_spread, oxygen_spread = (
        collect_measurements(bead_positions, symbols or [], cell, args.discard)
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for axis, values, title, xlabel in (
        (
            axes[0],
            bond_distances,
            "Intramolecular O-H bond",
            "O-H distance (angstrom)",
        ),
        (
            axes[1],
            neighbor_distances,
            "Nearest intermolecular O-H",
            "Nearest O-H distance (angstrom)",
        ),
    ):
        x, density = smooth_density(values)
        axis.plot(x, density, linewidth=2)
        axis.set(title=title, xlabel=xlabel, ylabel="Probability density (1/angstrom)")
        axis.grid(alpha=0.25)
    fig.savefig(args.output_dir / "distance_distributions.png", dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6, 4), constrained_layout=True)
    for element, values, color in (
        ("H", hydrogen_spread, "tab:blue"),
        ("O", oxygen_spread, "tab:red"),
    ):
        x, density = smooth_density(values)
        axis.plot(x, density, linewidth=2, color=color, label=element)
    axis.set(
        title="Ring-polymer bead spread from centroid",
        xlabel="Bead-to-centroid distance (angstrom)",
        ylabel="Probability density (1/angstrom)",
    )
    axis.grid(alpha=0.25)
    axis.legend(title="Element")
    fig.savefig(args.output_dir / "ring_polymer_spread.png", dpi=180)
    plt.close(fig)

    rows = [
        summarize("intramolecular_OH", bond_distances),
        summarize("nearest_inter_OH", neighbor_distances),
        summarize("hydrogen_bead_spread", hydrogen_spread),
        summarize("oxygen_bead_spread", oxygen_spread),
    ]
    with (args.output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Analyzed {len(trajectories)} beads and {len(bead_positions[0])} frames.")
    print(f"Results saved under {args.output_dir}.")


if __name__ == "__main__":
    main()
