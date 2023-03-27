# BioPandas
# License: BSD 3 clause
# Author goniochromatic: https://github.com/goniochromatic/
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
import os
import pandas as pd


TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), "data", "1t48_995.pdb")

p1t48 = PandasPdb()
p1t48.read_pdb(TESTDATA_1t48)


def test_accuracy():
    # Create test DataFrame with 3 atoms
    test_df = pd.DataFrame({"element_symbol": ["C", "O", "N", "S"],
                            "x_coord": [1.0, 2.0, 3.0, 4.0],
                            "y_coord": [5.0, 6.0, 7.0, 8.0],
                            "z_coord": [9.0, 10.0, 11.0, 12.0]})

    coords = test_df[["x_coord", "y_coord", "z_coord"]].to_numpy()
    masses = [12.0107, 15.9994, 14.0067, 32.065]
    total_mass = sum(masses)

    weighted_coords = [(m*x, m*y, m*z) for (x, y, z), m in zip(coords, masses)]
    weighted_deviation = sum(m * (x**2 + y**2 + z**2) for (x, y, z), m in zip(coords, masses))
    mean_weighted_coords = [sum(coords) / total_mass for coords in zip(*weighted_coords)]
    mean_weighted_deviation = sum(coord**2 for coord in mean_weighted_coords)

    # rounding needs to be changed to pass test if final rounding in gyradius will be changed
    expected_rg = round((weighted_deviation / total_mass - mean_weighted_deviation)**0.5, 4)
    rg = PandasPdb.gyradius(test_df)
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


def test_equal():
    rg = PandasPdb.gyradius(p1t48.df['ATOM'])
    expected_rg = 18.1508
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"
