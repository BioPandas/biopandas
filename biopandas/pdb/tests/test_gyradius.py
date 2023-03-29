# BioPandas
# License: BSD 3 clause
# Author goniochromatic: https://github.com/goniochromatic/
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
import os
import pandas as pd
from nose.tools import raises
import warnings

TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), "data", "1t48_995.pdb")

p1t48 = PandasPdb()
p1t48.read_pdb(TESTDATA_1t48)


def test_default_record():
    rg = p1t48.gyradius()
    expected_rg = 18.1508
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


def test_atom():
    rg = p1t48.gyradius(("ATOM",))
    expected_rg = 18.1508
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


def test_hetatm():
    rg = p1t48.gyradius(("HETATM",))
    expected_rg = 21.0143
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


def test_atom_and_hetatm():
    rg = p1t48.gyradius(("ATOM", "HETATM"))
    expected_rg = 19.0307
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


@raises(KeyError)
def test_wrong_record_name():
    p1t48.gyradius(("Wrong",))


@raises(TypeError)
def test_wrong_arg_type():
    p1t48.gyradius(5)


def test_decimals():
    rg = p1t48.gyradius(decimals=5)
    expected_rg = 18.15084
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


def test_negative_decimals():
    rg = p1t48.gyradius(decimals=-1)
    expected_rg = 20.0
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"


@raises(TypeError)
def test_wrong_decimals_arg():
    p1t48.gyradius(decimals='five')


def test_both_args():
    rg = p1t48.gyradius(("HETATM",), 5)
    expected_rg = 21.01426
    assert rg == expected_rg, f"Expected {expected_rg}, got {rg} instead"
