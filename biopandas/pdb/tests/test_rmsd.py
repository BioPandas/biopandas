# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
import os
from nose.tools import raises


TESTDATA_1t48 = os.path.join(os.path.dirname(__file__), 'data',
                                                        '1t48_995.pdb')
TESTDATA_1t49 = os.path.join(os.path.dirname(__file__), 'data',
                                                        '1t49_995.pdb')
TESTDATA_lig1 = os.path.join(os.path.dirname(__file__), 'data',
                                                        'lig_conf_1.pdb')
TESTDATA_lig2 = os.path.join(os.path.dirname(__file__), 'data',
                                                        'lig_conf_2.pdb')

p1t48 = PandasPdb()
p1t48.read_pdb(TESTDATA_1t48)
p1t49 = PandasPdb()
p1t49.read_pdb(TESTDATA_1t49)

pl1 = PandasPdb()
pl1.read_pdb(TESTDATA_lig1)
pl2 = PandasPdb()
pl2.read_pdb(TESTDATA_lig2)


def test_equal():
    r = PandasPdb.rmsd(p1t48.df['ATOM'], p1t48.df['ATOM'], s=None)
    assert r == 0.000, r


@raises(AttributeError)
def test_wrong_arg():
    PandasPdb.rmsd(p1t48.df['ATOM'].loc[1:, :], p1t48.df['ATOM'], s='bla')


@raises(AttributeError)
def test_incompatible():
    PandasPdb.rmsd(p1t48.df['ATOM'].loc[1:, :], p1t48.df['ATOM'], s=None)


@raises(AttributeError)
def test_invalid_query():
    PandasPdb.rmsd(p1t48.df['ATOM'].loc[1:, :], p1t48.df['ATOM'], s='bla')


def test_protein():
    r = PandasPdb.rmsd(p1t48.df['ATOM'], p1t49.df['ATOM'],
                       s='c-alpha', invert=False)
    assert r == 0.4785, r


def test_ligand():
    r = PandasPdb.rmsd(pl1.df['HETATM'], pl2.df['HETATM'],
                       s='hydrogen', invert=True)
    assert r == 1.9959, r


def test_ligand_default():
    r = PandasPdb.rmsd(pl1.df['HETATM'], pl2.df['HETATM'],
                       s=None)
    assert r == 2.6444, r
