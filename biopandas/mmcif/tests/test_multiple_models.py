# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# Author: Arian Jamasb <arian@jamasb.io>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas
import os

from biopandas.mmcif import PandasMmcif

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "2jyf.cif.gz")

def test_label_models():
    df = PandasMmcif().read_mmcif(TESTDATA_FILENAME)
    df.label_models()
    assert "model_id" in df.df["ATOM"].columns
    
def test_get_model():
    df = PandasMmcif().read_mmcif(TESTDATA_FILENAME)
    MODEL_INDEX = 1
    new_df = df.get_model(MODEL_INDEX)
    print(df)
    assert new_df.df["ATOM"]["pdbx_PDB_model_num"].all() == MODEL_INDEX


def test_get_models():
    df = PandasMmcif().read_mmcif(TESTDATA_FILENAME)
    MODEL_INDICES = [1, 3, 5]

    new_df = df.get_models(MODEL_INDICES)
    assert new_df.df["ATOM"]["pdbx_PDB_model_num"].all() in MODEL_INDICES