"""Generate the fixed-cell q-TIP4P/F LAMMPS structure from the PDB coordinates."""

from __future__ import annotations

from pathlib import Path

BOHR_PER_ANGSTROM = 1.889726124565
O_CHARGE = -1.1128
H_CHARGE = 0.5564


def read_water_box(path: Path) -> tuple[float, list[tuple[str, tuple[float, float, float]]]]:
    cell_angstrom: tuple[float, float, float] | None = None
    atoms: list[tuple[str, tuple[float, float, float]]] = []

    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("CRYST1"):
            a, b, c = float(line[6:15]), float(line[15:24]), float(line[24:33])
            alpha, beta, gamma = float(line[33:40]), float(line[40:47]), float(line[47:54])
            if any(abs(angle - 90.0) > 1e-5 for angle in (alpha, beta, gamma)):
                raise ValueError("The input PDB cell must be orthorhombic.")
            if max(a, b, c) - min(a, b, c) > 1e-5:
                raise ValueError("The input PDB cell must be cubic.")
            cell_angstrom = (a, b, c)
        elif line.startswith(("ATOM  ", "HETATM")):
            symbol = line[76:78].strip().upper()
            position = tuple(float(line[start : start + 8]) for start in (30, 38, 46))
            atoms.append((symbol, position))

    if cell_angstrom is None:
        raise ValueError(f"No CRYST1 unit cell found in {path}.")
    if len(atoms) == 0 or len(atoms) % 3 != 0:
        raise ValueError("Expected a non-empty set of complete three-atom water molecules.")
    for offset in range(0, len(atoms), 3):
        if [atom[0] for atom in atoms[offset : offset + 3]] != ["O", "H", "H"]:
            raise ValueError("Water atoms must be listed as O, H, H for each molecule.")

    cell_bohr = cell_angstrom[0] * BOHR_PER_ANGSTROM
    return cell_bohr, atoms


def wrapped_position(
    position: tuple[float, float, float], cell: float
) -> tuple[tuple[float, float, float], tuple[int, int, int]]:
    images = tuple(int(value // cell) for value in position)
    wrapped = tuple(value - image * cell for value, image in zip(position, images))
    return wrapped, images


def write_lammps_data(pdb_path: Path, output_path: Path) -> None:
    cell, atoms = read_water_box(pdb_path)
    molecule_count = len(atoms) // 3
    records: list[tuple[int, int, int, float, tuple[float, float, float], tuple[int, int, int]]] = []

    for molecule_index in range(molecule_count):
        start = molecule_index * 3
        oxygen = tuple(value * BOHR_PER_ANGSTROM for value in atoms[start][1])
        oxygen_wrapped, oxygen_image = wrapped_position(oxygen, cell)
        records.append((molecule_index + 1, molecule_index + 1, 1, O_CHARGE, oxygen_wrapped, oxygen_image))

        for hydrogen_offset in (1, 2):
            hydrogen = tuple(
                value * BOHR_PER_ANGSTROM for value in atoms[start + hydrogen_offset][1]
            )
            displacement = tuple(
                h - o - round((h - o) / cell) * cell
                for h, o in zip(hydrogen, oxygen)
            )
            same_image_position = tuple(
                oxygen_wrapped[axis] + displacement[axis] for axis in range(3)
            )
            hydrogen_wrapped, hydrogen_image = wrapped_position(same_image_position, cell)
            atom_id = start + hydrogen_offset + 1
            records.append(
                (
                    atom_id,
                    molecule_index + 1,
                    2,
                    H_CHARGE,
                    hydrogen_wrapped,
                    hydrogen_image,
                )
            )

    lines = [
        "q-TIP4P/F water box generated from data/water_32.pdb",
        "",
        f"{len(atoms)} atoms",
        f"{2 * molecule_count} bonds",
        f"{molecule_count} angles",
        "",
        "2 atom types",
        "1 bond types",
        "1 angle types",
        "",
        f"0.0 {cell:.8f} xlo xhi",
        f"0.0 {cell:.8f} ylo yhi",
        f"0.0 {cell:.8f} zlo zhi",
        "",
        "Masses",
        "",
        "1 15.9994",
        "2 1.0080",
        "",
        "Bond Coeffs",
        "",
        "1 1.78 0.2708585 -0.327738785 0.231328959",
        "",
        "Angle Coeffs",
        "",
        "1 0.0700 107.400000",
        "",
        "Atoms # full",
        "",
    ]

    for atom_id, molecule_id, atom_type, charge, position, images in records:
        x, y, z = position
        nx, ny, nz = images
        lines.append(
            f"{atom_id} {molecule_id} {atom_type} {charge:.4f} "
            f"{x:.8f} {y:.8f} {z:.8f} {nx} {ny} {nz}"
        )

    lines.extend(("", "Bonds", ""))
    for molecule_index in range(molecule_count):
        oxygen_id = 3 * molecule_index + 1
        lines.append(f"{2 * molecule_index + 1} 1 {oxygen_id} {oxygen_id + 1}")
        lines.append(f"{2 * molecule_index + 2} 1 {oxygen_id} {oxygen_id + 2}")

    lines.extend(("", "Angles", ""))
    for molecule_index in range(molecule_count):
        oxygen_id = 3 * molecule_index + 1
        lines.append(f"{molecule_index + 1} 1 {oxygen_id + 1} {oxygen_id} {oxygen_id + 2}")

    output_path.write_text("\n".join(lines) + "\n", encoding="ascii")


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parent.parent
    write_lammps_data(
        project_dir / "data" / "water_32.pdb",
        project_dir / "data" / "water_32_data.lmp",
    )
