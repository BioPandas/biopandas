# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import numpy as np
from biopandas.pdb import PandasPdb
import os


def test_defaults():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), 'data',
                                                            '1t48_995.pdb')
    p1t48 = PandasPdb()
    p1t48.read_pdb(TESTDATA_1t48)
    expect_res = ['M', 'E', 'M', 'E', 'K', 'E', 'F', 'E', 'Q',
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

    transl = p1t48.amino3to1()
    expect_chain = ['A' for _ in range(transl.shape[0])]
    got_chain = list(transl['chain_id'].values)
    got_res = list(transl['residue_name'].values)

    assert expect_chain == got_chain
    assert expect_res == got_res


def test_sameindex():
    TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), 'data',
                                                            '1t48_995.pdb')
    p1t48 = PandasPdb()
    p1t48.read_pdb(TESTDATA_1t48)
    print(p1t48)
    p1t48.df['ATOM'].index = np.zeros(p1t48.df['ATOM'].shape[0], dtype=int)

    expect_res = ['M', 'E', 'M', 'E', 'K', 'E', 'F', 'E', 'Q',
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

    transl = p1t48.amino3to1()
    expect_chain = ['A' for _ in range(transl.shape[0])]
    got_chain = list(transl['chain_id'].values)
    got_res = list(transl['residue_name'].values)

    assert expect_chain == got_chain
    assert expect_res == got_res


def test_multichain():
    TESTDATA_5mtn = os.path.join(os.path.dirname(__file__),
                                 'data', '5mtn_multichain.pdb')
    mtn = PandasPdb()
    mtn.read_pdb(TESTDATA_5mtn)
    expect_res_a = ['S', 'L', 'E', 'P', 'E', 'P', 'W', 'F', 'F', 'K', 'N', 'L',
                    'S', 'R', 'K', 'D', 'A', 'E', 'R', 'Q', 'L', 'L', 'A', 'P',
                    'G', 'N', 'T', 'H', 'G', 'S', 'F', 'L', 'I', 'R', 'E', 'S',
                    'E', 'S', 'T', 'A', 'G', 'S', 'F', 'S', 'L', 'S', 'V', 'R',
                    'D', 'F', 'D', 'Q', 'G', 'E', 'V', 'V', 'K', 'H', 'Y', 'K',
                    'I', 'R', 'N', 'L', 'D', 'N', 'G', 'G', 'F', 'Y', 'I', 'S',
                    'P', 'R', 'I', 'T', 'F', 'P', 'G', 'L', 'H', 'E', 'L', 'V',
                    'R', 'H', 'Y', 'T']
    expect_res_b = ['S', 'V', 'S', 'S', 'V', 'P', 'T', 'K', 'L', 'E', 'V', 'V',
                    'A', 'A', 'T', 'P', 'T', 'S', 'L', 'L', 'I', 'S', 'W', 'D',
                    'A', 'P', 'A', 'V', 'T', 'V', 'V', 'Y', 'Y', 'L', 'I', 'T',
                    'Y', 'G', 'E', 'T', 'G', 'S', 'P', 'W', 'P', 'G', 'G', 'Q',
                    'A', 'F', 'E', 'V', 'P', 'G', 'S', 'K', 'S', 'T', 'A', 'T',
                    'I', 'S', 'G', 'L', 'K', 'P', 'G', 'V', 'D', 'Y', 'T', 'I',
                    'T', 'V', 'Y', 'A', 'H', 'R', 'S', 'S', 'Y', 'G', 'Y', 'S',
                    'E', 'N', 'P', 'I', 'S', 'I', 'N', 'Y', 'R', 'T']

    transl = mtn.amino3to1()

    expect_chain = ['A' for _ in range(88)] + ['B' for _ in range(94)]
    got_chain = list(transl['chain_id'].values)

    got_res_a = list(transl.loc[transl['chain_id'] == 'A',
                                'residue_name'].values)
    got_res_b = list(transl.loc[transl['chain_id'] == 'B',
                                'residue_name'].values)

    assert expect_chain == got_chain
    assert expect_res_a == got_res_a
    assert expect_res_b == got_res_b


def test_pdb_with_insertion_codes():

    PDB_2D7T_PATH = os.path.join(os.path.dirname(__file__),
                                 'data', '2d7t.pdb')

    ppdb = PandasPdb().read_pdb(PDB_2D7T_PATH)
    sequence = ppdb.amino3to1()
    assert "".join(sequence[50:60]['residue_name'].values) == 'INPKSGDTNY'
