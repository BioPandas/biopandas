# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
from biopandas.testutils import assert_raises
import os


TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'data', '3eiy.pdb')


def test_overwrite_df():
    data_path = os.path.join(os.path.dirname(__file__), 'data', '3eiy.pdb')
    pdb = PandasPdb().read_pdb(data_path)

    def overwrite():
        pdb.df = 'bla'

    expect = ('Please use `PandasPdb._df = ... ` instead\n'
              'of `PandasPdb.df = ... ` if you are sure that\n'
              'you want to overwrite the `df` attribute.')

    assert_raises(AttributeError,
                  expect,
                  overwrite)
