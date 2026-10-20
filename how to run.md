# How to run the water PIMD simulation

This guide walks through setup, running the simulation, and plotting its output.
The commands below use **WSL2 (Ubuntu)** because the LAMMPS/i-PI connection
uses a Unix-domain socket. Run them from the project folder inside WSL, not
from native PowerShell.

## 1. Open the project in WSL

If the project is on your Windows `C:` drive, open Ubuntu/WSL and change to it:

```bash
cd "/mnt/c/Users/<your-windows-username>/GitHub Projects/Effects_Of_Water_Quantum_Hydrogen_effects"
```

If you cloned the project into your WSL home directory instead, `cd` to that
clone's directory.

## 2. Install the software

Install Miniforge or another Conda distribution in WSL if you do not already
have one. Then create and activate the project environment:

```bash
conda env create -f environment.yml
conda activate water-pimd
```

Confirm that i-PI, LAMMPS, and the required LAMMPS styles are available:

```bash
command -v i-pi
command -v lmp
lmp -h
```

This simulation needs a LAMMPS build with the `MISC` package (`fix ipi`),
`CLASS2` (`bond_style class2`), `KSPACE` (PPPM/TIP4P), and `MOLECULE`
support. Some Conda builds may not include every optional package. If `lmp -h`
does not list the needed styles, install/build a LAMMPS version that includes
those packages before continuing.

## 3. Generate the LAMMPS water-box data

From the project root, run:

```bash
python data/prepare_lammps_data.py
```

This reads the periodic 32-water structure in `data/water_32.pdb` and writes
`data/water_32_data.lmp`, which contains the matching cell, charges, masses,
and q-TIP4P/F bond/angle topology.

## 4. Start the simulation

Run:

```bash
bash run.sh
```

The script starts the i-PI server, then two LAMMPS force clients. They connect
over the local Unix socket named `h2o-lammps`. The current settings in
`data/input.xml` are:

- 32 path-integral beads
- 298 K
- 0.25 fs timestep
- 4,000 steps (1 ps total; a short smoke run only)

Progress is written to `simulation.log` and `lammps-1.log` /
`lammps-2.log`. The simulation outputs have the `simulation` prefix. A
successful run should create `simulation.out` and bead trajectories named
`simulation.pos_0.xyz` through `simulation.pos_31.xyz`.

If a previous run left an i-PI process or socket behind, stop that process
before starting another run. Do not run two copies in the same directory at
the same time.

## 5. Plot and summarize the results

After the simulation completes, run:

```bash
python analyze.py --discard 20
```

`--discard 20` skips the first 20 saved frames as an initial transient. Adjust
the count based on equilibration; for a short smoke run, there may not be
enough production frames for meaningful averages.

The analysis writes:

- `analysis/distance_distributions.png` — intramolecular O-H bond distances
  and nearest intermolecular O-H distances.
- `analysis/ring_polymer_spread.png` — hydrogen and oxygen bead distances
  from their ring-polymer centroids.
- `analysis/summary.csv` — means, medians, and 95th percentiles.

The analysis expects the atom order in each XYZ frame to be O-H-H for every
water molecule, matching the provided input structure.

## 6. Prepare a longer, converged calculation

The default 4,000-step calculation is for verifying that the workflow runs; it
is not enough to establish converged water properties. Before interpreting
results:

1. Equilibrate the system, then gather substantially longer production
   trajectories (often at least hundreds of picoseconds, depending on the
   property and sampling uncertainty).
2. Compare runs with 32 and 64 beads, increasing `nbeads` in `data/input.xml`.
3. Check timestep stability, independent seeds, and sampling uncertainty.
4. Validate the chosen q-TIP4P/F model against appropriate reference data for
   the property you are studying.

Increase `<total_steps>` in `data/input.xml` to extend a run. For example,
400,000 steps at 0.25 fs correspond to 100 ps. Larger bead counts and long
trajectories increase computational cost.

PIMD samples equilibrium nuclear quantum effects. Its thermostat-driven
trajectory time is not a direct physical real-time vibration timescale.
Broadened distributions and visible hydrogen-bead spread are useful
diagnostics, but do not alone prove that a particular observable is converged
or quantitatively accurate.
