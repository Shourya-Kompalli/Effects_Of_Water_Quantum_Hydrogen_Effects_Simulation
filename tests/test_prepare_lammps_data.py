import tempfile
import unittest
from pathlib import Path

from data.prepare_lammps_data import write_lammps_data


PROJECT_DIR = Path(__file__).resolve().parents[1]


class PrepareLammpsDataTests(unittest.TestCase):
    def test_writes_periodic_neutral_32_water_topology(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_path = Path(temporary_directory) / "water.lmp"
            write_lammps_data(PROJECT_DIR / "data" / "water_32.pdb", output_path)
            content = output_path.read_text(encoding="ascii")

        self.assertIn("96 atoms", content)
        self.assertIn("64 bonds", content)
        self.assertIn("32 angles", content)
        self.assertIn("0.0 19.38859", content)

        atom_section = content.split("Atoms # full\n\n", 1)[1].split("\n\nBonds", 1)[0]
        atom_records = [line.split() for line in atom_section.splitlines()]
        self.assertEqual(len(atom_records), 96)
        self.assertAlmostEqual(sum(float(record[3]) for record in atom_records), 0.0)

        cell = 19.38859
        for record in atom_records:
            for coordinate in record[4:7]:
                self.assertGreaterEqual(float(coordinate), 0.0)
                self.assertLess(float(coordinate), cell)


if __name__ == "__main__":
    unittest.main()
