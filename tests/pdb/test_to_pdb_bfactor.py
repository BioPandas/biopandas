"""Regression test for #109: b_factor must fit the 6-char PDB field.

ATOM b_factor is fixed at columns 61-66 (6 chars). Values that round up
across an integer-digit boundary (e.g. 99.996 -> 100.00) must not spill
into the segment_id/element columns (cols 73-78).
"""
import os
import tempfile
import unittest

import pandas as pd

from biopandas.pdb import PandasPdb


def _write_pdb(b_factor):
    df = pd.DataFrame([{
        "record_name": "ATOM",
        "atom_number": 1,
        "blank_1": " ",
        "atom_name": "CA",
        "alt_loc": " ",
        "residue_name": "ASP",
        "blank_2": " ",
        "chain_id": "A",
        "residue_number": 100,
        "insertion": " ",
        "blank_3": "   ",
        "x_coord": 219.123,
        "y_coord": 233.404,
        "z_coord": 332.880,
        "occupancy": 1.00,
        "b_factor": b_factor,
        "blank_4": "      ",
        "segment_id": "",
        "element_symbol": "C",
        "charge": 0,
        "line_idx": 1,
    }])
    pp = PandasPdb()
    pp._df = {"ATOM": df}
    path = tempfile.mktemp(suffix=".pdb")
    pp.to_pdb(path)
    with open(path) as f:
        line = f.readline().rstrip("\n")
    os.remove(path)
    return line


class TestBFactorFieldWidth(unittest.TestCase):
    def test_bfactor_fits_six_char_field(self):
        # all of these fit in the 6-char b_factor field (cols 61-66);
        # 99.996 is the rounding-up edge (formats as 100.00)
        for b in (97.39, 99.99, 99.996, 100.32, 999.99):
            line = _write_pdb(b)
            self.assertEqual(
                len(line), 80,
                "b_factor=%s produced a %d-char ATOM line" % (b, len(line)),
            )
            self.assertEqual(
                line[76:78], " C",
                "element column shifted for b_factor=%s: %r" % (b, line),
            )


if __name__ == "__main__":
    unittest.main()
