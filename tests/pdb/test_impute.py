# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas


import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import tests.pdb.data
from biopandas.pdb import PandasPdb

TEST_DATA = pkg_resources.files(tests.pdb.data)


TESTDATA_FILENAME = str(TEST_DATA.joinpath("3eiy_stripped_no_ele.pdb"))

ppdb = PandasPdb()
ppdb.read_pdb(TESTDATA_FILENAME)


def test_impute_hetatm():
    new = ppdb.impute_element(records=["HETATM"])
    assert new["HETATM"]["element_symbol"][1] == "N"
    assert new["HETATM"]["element_symbol"][10] == "O"
    assert new["ATOM"]["element_symbol"][1] == ""
    assert new["ATOM"]["element_symbol"][10] == ""


def test_impute_atom():
    new = ppdb.impute_element(records=["ATOM"])
    assert new["ATOM"]["element_symbol"][1] == "C"
    assert new["ATOM"]["element_symbol"][10] == "C"
    assert new["HETATM"]["element_symbol"][1] == ""
    assert new["HETATM"]["element_symbol"][10] == ""


def test_impute_default():
    assert ppdb.df["ATOM"]["element_symbol"][1] == ""
    assert ppdb.df["HETATM"]["element_symbol"][1] == ""
    new = ppdb.impute_element()
    assert new["ATOM"]["element_symbol"][1] == "C"
    assert new["ATOM"]["element_symbol"][10] == "C"
    assert new["HETATM"]["element_symbol"][1] == "N"
    assert new["HETATM"]["element_symbol"][10] == "O"
