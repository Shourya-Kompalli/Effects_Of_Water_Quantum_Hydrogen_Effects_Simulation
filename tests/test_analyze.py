import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from analyze import collect_measurements, smooth_density

PROJECT_DIR = Path(__file__).resolve().parents[1]


class AnalyzeTests(unittest.TestCase):
    def test_collects_distances_and_element_specific_bead_spread(self) -> None:
        symbols = ["O", "H", "H", "O", "H", "H"]
        water_pair = np.asarray(
            [
                [1.0, 1.0, 1.0],
                [1.96, 1.0, 1.0],
                [1.0, 1.96, 1.0],
                [4.0, 1.0, 1.0],
                [4.96, 1.0, 1.0],
                [4.0, 1.96, 1.0],
            ]
        )
        positions = np.tile(water_pair, (2, 2, 1, 1))
        positions[1, :, 1, 0] += 0.2

        bonds, neighbors, hydrogen_spread, oxygen_spread = collect_measurements(
            [positions[0], positions[1]], symbols, cell=10.0, discard=0
        )

        self.assertEqual(bonds.size, 16)
        self.assertEqual(neighbors.size, 16)
        self.assertTrue(np.all(np.isfinite(neighbors)))
        self.assertGreater(np.std(bonds), 0.0)
        self.assertGreater(np.max(hydrogen_spread), 0.0)
        self.assertTrue(np.allclose(oxygen_spread, 0.0))

    def test_smoothed_density_is_normalized(self) -> None:
        values = np.linspace(0.8, 1.2, 100)
        x, density = smooth_density(values)

        self.assertEqual(x.size, 120)
        self.assertAlmostEqual(float(np.sum(density) * (x[1] - x[0])), 1.0)

    def test_cli_writes_distance_plots_and_summary(self) -> None:
        symbols = ["O", "H", "H", "O", "H", "H"]
        water_pair = np.asarray(
            [
                [1.0, 1.0, 1.0],
                [1.96, 1.0, 1.0],
                [1.0, 1.96, 1.0],
                [4.0, 1.0, 1.0],
                [4.96, 1.0, 1.0],
                [4.0, 1.96, 1.0],
            ]
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            for bead in range(32):
                trajectory = temporary_path / f"synthetic.pos_{bead}.xyz"
                with trajectory.open("w", encoding="utf-8") as output:
                    for frame in range(3):
                        positions = water_pair.copy()
                        hydrogen_shift = 0.1 * np.sin(2.0 * np.pi * (bead + frame) / 32)
                        positions[1, 0] += hydrogen_shift
                        positions[4, 1] += hydrogen_shift
                        output.write(f"{len(symbols)}\nframe {frame}\n")
                        for symbol, position in zip(symbols, positions):
                            output.write(
                                f"{symbol} {position[0]} {position[1]} {position[2]}\n"
                            )

            output_dir = temporary_path / "analysis"
            subprocess.run(
                [
                    sys.executable,
                    str(PROJECT_DIR / "analyze.py"),
                    "--prefix",
                    "synthetic",
                    "--structure",
                    str(PROJECT_DIR / "data" / "water_32.pdb"),
                    "--output-dir",
                    str(output_dir),
                    "--discard",
                    "0",
                ],
                cwd=temporary_path,
                check=True,
                capture_output=True,
                text=True,
            )

            self.assertGreater((output_dir / "distance_distributions.png").stat().st_size, 0)
            self.assertGreater((output_dir / "ring_polymer_spread.png").stat().st_size, 0)
            with (output_dir / "summary.csv").open(encoding="utf-8") as summary_file:
                rows = list(csv.DictReader(summary_file))
            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[2]["measurement"], "hydrogen_bead_spread")


if __name__ == "__main__":
    unittest.main()
