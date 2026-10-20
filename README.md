# Quantum water PIMD starter

This project is a runnable starting point for liquid-water path-integral
molecular dynamics (PIMD). It uses i-PI for the ring-polymer dynamics and
LAMMPS for the flexible q-TIP4P/F water potential. The initial structure is a
periodic 32-molecule water box.

**DONE BY VIBE CODING WITH COPIOLOT**

## Why this setup

- Each nucleus is represented by 32 imaginary-time beads at 298 K. Increase
  this to 64 and compare observables to check bead convergence.
- q-TIP4P/F is a flexible, anharmonic four-site water model. Its electrostatic
  charge site is displaced from oxygen; hydrogen-bond effects arise from the
  intermolecular force field, not from a separate "quantum tug" term.
- The box and all interaction parameters are defined for a fixed periodic cell
  on the LAMMPS side. i-PI controls the PIMD sampling.
- The default run is a short smoke run, not production data. The analysis
  script reports intramolecular O-H distances and each element's bead-to-
  centroid spread.

LAMMPS `fix ipi` uses Unix-domain sockets, so run this setup in Linux or WSL2,
not native Windows. The published LAMMPS/i-PI q-TIP4P/F example that informs
the force-field coefficients and client pattern is in the
[atomistic-cookbook path-integrals example](https://github.com/lab-cosmo/atomistic-cookbook/tree/fff8b352d5ff6028a26cbb60f1abe7b07b8a968d/examples/path-integrals).
The water coordinates are from that example; `data/prepare_lammps_data.py`
generates a LAMMPS data file with the matching periodic cell and topology.

## Install (WSL2)

From the project directory in WSL:

```bash
conda env create -f environment.yml
conda activate water-pimd
```

For complete setup, run, and analysis instructions, see [how to run.md](how%20to%20run.md).

Check that the installed LAMMPS build includes the `ipi` fix and the `class2`,
`kspace`, and `molecule` packages:

```bash
lmp -h
```

If `fix ipi` or `bond_style class2` is unavailable, use a LAMMPS build compiled
with the MISC, CLASS2, KSPACE, and MOLECULE packages before running.

## Run

```bash
bash run.sh
```

The script starts the i-PI server and two LAMMPS force clients. Output files
use the `simulation` prefix. For a quick health check, confirm the i-PI log
reaches the configured step count and that `simulation.pos_0.xyz` through
`simulation.pos_31.xyz` are present.

The default is 4,000 steps at 0.25 fs (1 ps), which is only enough to exercise
the workflow. For production, equilibrate first, then use at least hundreds of
picoseconds of sampling and check independent runs and convergence with bead
count, timestep, system size, and sampling duration. For example, 400,000
steps at 0.25 fs is 100 ps; edit `<total_steps>` in `data/input.xml` and
consider a 64-bead comparison. These runs can be computationally expensive.

## Analyze

After a run:

```bash
python analyze.py --discard 20
```

This writes:

- `analysis/distance_distributions.png`: intramolecular O-H bond distances and
  each hydrogen's nearest intermolecular oxygen distance.
- `analysis/ring_polymer_spread.png`: distributions of bead distances from
  each atom's ring-polymer centroid, separated by element.
- `analysis/summary.csv`: means and percentiles for these measurements.

Use enough equilibrated frames to obtain smooth, stable distributions. A broad
curve alone does not establish physical accuracy: compare the 32- and 64-bead
results and validate against appropriate experimental or higher-level
reference data. PIMD trajectories sample equilibrium nuclear quantum effects;
their thermostat-driven time axis is not a direct measure of real-time
vibrational frequencies.

## Important model boundary

This is a classical empirical potential with quantum nuclei, not an electronic
structure calculation. q-TIP4P/F is a useful water model for demonstrating
nuclear quantum effects, but quantitative conclusions still require
convergence tests and validation for the specific property being studied.
