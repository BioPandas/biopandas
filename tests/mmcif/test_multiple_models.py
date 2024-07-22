# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# Author: Arian Jamasb <arian@jamasb.io>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas
import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import tests.mmcif.data
from biopandas.mmcif import PandasMmcif

TEST_DATA = pkg_resources.files(tests.mmcif.data)

TESTDATA_FILENAME = str(TEST_DATA.joinpath("2jyf.cif.gz"))


def test_label_models():
    biopandas_structure = PandasMmcif().read_mmcif(TESTDATA_FILENAME)
    biopandas_structure.label_models()
    assert "model_id" in biopandas_structure.df["ATOM"].columns


def test_get_model():
    biopandas_structure = PandasMmcif().read_mmcif(TESTDATA_FILENAME)
    MODEL_INDEX = 1
    new_biopandas_structure = biopandas_structure.get_model(MODEL_INDEX)
    assert (
        new_biopandas_structure.df["ATOM"]["pdbx_PDB_model_num"].all()
        == MODEL_INDEX
    )


def test_get_models():
    biopandas_structure = PandasMmcif().read_mmcif(TESTDATA_FILENAME)
    MODEL_INDICES = [1, 3, 5]

    new_biopandas_structure = biopandas_structure.get_models(MODEL_INDICES)
    assert (
        new_biopandas_structure.df["ATOM"]["pdbx_PDB_model_num"].all()
        in MODEL_INDICES
    )
