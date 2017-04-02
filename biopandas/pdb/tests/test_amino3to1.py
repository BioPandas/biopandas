# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
import os


def test_defaults():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), 'data',
                                                            '1t48_995.pdb')
    p1t48 = PandasPdb()
    p1t48.read_pdb(TESTDATA_1t48)
    expect = ['M', 'E', 'M', 'E', 'K', 'E', 'F', 'E', 'Q',
              'I', 'D', 'K', 'S', 'G', 'S', 'W', 'A', 'A',
              'I', 'Y', 'Q', 'D', 'I', 'R', 'H', 'E', 'A',
              'S', 'D', 'F', 'P', 'C', 'R', 'V', 'A', 'K',
              'L', 'P', 'K', 'N', 'K', 'N', 'R', 'N', 'R',
              'Y', 'R', 'D', 'V', 'S', 'P', 'F', 'D', 'H',
              'S', 'R', 'I', 'K', 'L', 'H', 'Q', 'E', 'D',
              'N', 'D', 'Y', 'I', 'N', 'A', 'S', 'L', 'I',
              'K', 'M', 'E', 'E', 'A', 'Q', 'R', 'S', 'Y',
              'I', 'L', 'T', 'Q', 'G', 'P', 'L', 'P', 'N',
              'T', 'C', 'G', 'H', 'F', 'W', 'E', 'M', 'V',
              'W', 'E', 'Q', 'K', 'S', 'R', 'G', 'V', 'V',
              'M', 'L', 'N', 'R', 'V', 'M', 'E', 'K', 'G',
              'S', 'L', 'K']
    assert expect == list(p1t48.amino3to1().values)
