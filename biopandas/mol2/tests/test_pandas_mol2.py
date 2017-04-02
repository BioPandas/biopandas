# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import os
from biopandas.mol2 import PandasMOL2
from biopandas.mol2.mol2_io import split_multimol2

this_dir = os.path.dirname(os.path.realpath(__file__))


def test_read_mol2():

    data_path_1 = os.path.join(this_dir, 'data', '40_mol2_files.mol2')
    data_path_2 = os.path.join(this_dir, 'data', '40_mol2_files.mol2.gz')

    for data_path in (data_path_1, data_path_2):
        pdmol = PandasMOL2().read_mol2(data_path)
        assert pdmol.df.shape == (65, 9)
        assert pdmol.code == 'ZINC38611810'

        expect = ['atom_id', 'atom_name', 'x', 'y', 'z',
                  'atom_type', 'subst_id', 'subst_name', 'charge']
        assert expect == list(pdmol.df.columns)
        assert len(pdmol.mol2_text) == 6469


def test_read_mol2_from_list():

    data_path = os.path.join(this_dir, 'data', '40_mol2_files.mol2')
    mol2 = next(split_multimol2(data_path))

    pdmol = PandasMOL2().read_mol2_from_list(mol2_lines=mol2[1],
                                             mol2_code=mol2[0])
    assert pdmol.df.shape == (65, 9)
    assert pdmol.code == 'ZINC38611810'


def test_rmsd():
    data_path_1 = os.path.join(this_dir, 'data', '1b5e_1.mol2')
    data_path_2 = os.path.join(this_dir, 'data', '1b5e_2.mol2')

    pdmol_1 = PandasMOL2().read_mol2(data_path_1)
    pdmol_2 = PandasMOL2().read_mol2(data_path_2)

    assert pdmol_1.rmsd(pdmol_1.df, pdmol_2.df, heavy_only=False) == 1.5523
    assert pdmol_1.rmsd(pdmol_1.df, pdmol_2.df) == 1.1609


def test_distance():
    data_path = os.path.join(this_dir, 'data', '1b5e_1.mol2')

    pdmol = PandasMOL2().read_mol2(data_path)
    assert round(pdmol.distance().values[0], 3) == 31.185
