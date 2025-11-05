# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import os
import sys

import numpy as np
from nose.tools import assert_raises

from biopandas.align import TMAlign, Align
from biopandas.pdb import PandasPdb
from biopandas.stack.stack import PandasPdbStack

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")
TESTDATA_FILENAME2 = os.path.join(
    os.path.dirname(__file__), "data", "4eiy_anisouchunk.pdb"
)
TESTDATA_FILENAME3 = os.path.join(
    os.path.dirname(__file__), "data", "2jyf.pdb"
)
TESTDATA_FILENAME4 = os.path.join(
    os.path.dirname(__file__), "data", "2d7t.pdb"
)

OUTFILE = os.path.join(os.path.dirname(__file__), "data", "tmp.pdb")

def test_init():
    tmalign = TMAlign()

    # if windows tmalign.exe is present
    if sys.platform == 'win32':
        assert 'USalign.exe' in tmalign.tmalign_path
    else:
        assert 'USalign' in tmalign.tmalign_path

def test_init_with_path_ok():
    if sys.platform == 'win32':
        tmalign = TMAlign(os.path.join(os.path.dirname(__file__), "../../biopandas/align/USalign.exe"))
    else:
        tmalign = TMAlign(os.path.join(os.path.dirname(__file__), "../../biopandas/align/USalign"))

    assert 'USalign' in tmalign.tmalign_path

def test_init_with_path_not_ok():
    assert_raises(FileNotFoundError, TMAlign, 'fake_Usalign')

def test_init_with_path_none():
    if os.path.exists(os.path.join(os.path.dirname(__file__), "../../biopandas/align/USalign.exe")):
        tmalign = TMAlign()
        assert 'USalign' in tmalign.tmalign_path
    elif os.path.exists(os.path.join(os.path.dirname(__file__), "../../biopandas/align/USalign")):
        tmalign = TMAlign()
        assert 'USalign' in tmalign.tmalign_path
    else:
        assert_raises(ValueError, TMAlign, None)

def test_parse_tmalign_rotation_matrix():
    tmalign = TMAlign()
    matrix, translation = tmalign.parse_tmalign_rotation_matrix(
        os.path.join(os.path.dirname(__file__), "data", "tmalign_output.txt")
    )
    print("\n", matrix)
    assert matrix.shape == (3, 3)
    assert translation.shape == (3,)
    assert matrix[0, 0] == 0.9999934600
    assert matrix[0, 1] == 0.0029980568

def test_parse_tmalign_rotation_matrix_file_not_found():
    tmalign = TMAlign()
    assert_raises(FileNotFoundError, tmalign.parse_tmalign_rotation_matrix, 'file_not_found.txt')

def test_run_tmalign():
    tmalign = TMAlign()
    matrix_file_path, tm_score = tmalign.run_tmalign(TESTDATA_FILENAME, TESTDATA_FILENAME)

    assert matrix_file_path.endswith('.txt')
    assert tm_score == 1

def test_run_tmalign_file_not_found():
    tmalign = TMAlign()
    assert_raises(FileNotFoundError, tmalign.run_tmalign, TESTDATA_FILENAME, 'file_not_found.txt')
    assert_raises(FileNotFoundError, tmalign.run_tmalign, 'file_not_found.txt', TESTDATA_FILENAME)

def test_run_tmalign_structure_too_short():
    tmalign = TMAlign()
    assert_raises(ValueError, tmalign.run_tmalign, TESTDATA_FILENAME, TESTDATA_FILENAME2)

def test_process_structure_for_tmalign_perfect():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb2 = PandasPdb()
    ppdb2.read_pdb(TESTDATA_FILENAME)

    # apply random rotation
    random_rotation = np.array([[0.5, -0.866, 0],
                                [0.866, 0.5, 0],
                                [0, 0, 1]])

    random_translation = np.array([1, 2, 3])

    align = Align()
    coords_mobile = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == 'A'][['x_coord', 'y_coord', 'z_coord']].values
    transformed_coords = align.transform(coords_mobile, random_rotation, random_translation)
    ppdb2.df['ATOM'][['x_coord', 'y_coord', 'z_coord']] = transformed_coords

    target_file = tmalign.write_pdb_to_temp_file(ppdb)
    transformed_mobile, tm_score = tmalign.process_structure_for_tmalign(target_file, ppdb2, 'A')
    coords_transformed = transformed_mobile.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].values

    assert tm_score == 1
    assert np.allclose(coords_transformed, coords_mobile, atol=1e-3)

def test_process_structure_for_tmalign_protein_rna():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb2 = PandasPdb()
    ppdb2.read_pdb(TESTDATA_FILENAME3)

    target_file = tmalign.write_pdb_to_temp_file(ppdb)
    transformed_mobile, tm_score = tmalign.process_structure_for_tmalign(target_file, ppdb2, 'B')

    assert tm_score == 0.16383

def test_tmalign_to_stack_chains_exist():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME3, TESTDATA_FILENAME4])

    transformed_structures, tm_scores = tmalign.tmalign_to(ppdb, ppdb_stack,
                                                           'A',
                                                           {'2jyf': 'A', '2d7t': 'H', '3eiy': 'A'})
    assert tm_scores['2jyf'] == 0.16382
    assert tm_scores['2d7t'] == 0.32871
    assert tm_scores['3eiy'] == 1

def test_tmalign_to_stack_chains_missing():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME3, TESTDATA_FILENAME4])

    assert_raises(ValueError, tmalign.tmalign_to, ppdb, ppdb_stack, 'A', {'2jyf': 'A', '2d7t': 'H'})

def test_tmalign_to_stack_chains_not_exist():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs([TESTDATA_FILENAME3, TESTDATA_FILENAME4])

    assert_raises(ValueError, tmalign.tmalign_to, ppdb, ppdb_stack, 'A', {'2jyf': 'A', '2d7t': 'A'})

def test_tmalign_to_one():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb2 = PandasPdb()
    ppdb2.read_pdb(TESTDATA_FILENAME3)

    transformed_structures, tm_score = tmalign.tmalign_to(ppdb, ppdb2, 'A', 'A')
    assert tm_score == 0.16382

def test_tmalign_chain_missing():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb2 = PandasPdb()
    ppdb2.read_pdb(TESTDATA_FILENAME3)

    assert_raises(ValueError, tmalign.tmalign_to, ppdb, ppdb2, 'A', 'H')
    assert_raises(ValueError, tmalign.tmalign_to, ppdb, ppdb2, 'B', 'A')

def test_tmalign_to_wrong_input():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb2 = PandasPdb()
    ppdb2.read_pdb(TESTDATA_FILENAME3)

    assert_raises(ValueError, tmalign.tmalign_to, ppdb, 'random_string', 'A', 'A')

def test_tmalign_in_stack():
    tmalign = TMAlign()

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs(['1ycr', '2d7t', '3eiy', '4LWV'])

    chains = {'1ycr': 'A', '2d7t': 'H', '3eiy': 'A', '4LWV': 'A'}

    target_pdb_id, transformed_structures, tm_scores = tmalign.tmalign_in_stack(ppdb_stack, chains)

    assert tm_scores == {'2d7t': 0.27812, '3eiy': 0.23483, '4LWV': 0.92383}
    assert target_pdb_id == '1ycr'

def test_tmalign_in_stack_define_target():
    tmalign = TMAlign()

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs(['1ycr', '2d7t', '3eiy', '4LWV'])

    chains = {'1ycr': 'A', '2d7t': 'H', '3eiy': 'A', '4LWV': 'A'}

    target_pdb_id, transformed_structures, tm_scores = tmalign.tmalign_in_stack(ppdb_stack, chains, target='2d7t')

    assert tm_scores == {'1ycr': 0.33733, '3eiy': 0.24401, '4LWV': 0.30331}
    assert target_pdb_id == '2d7t'

def test_tmalign_in_stack_define_target_not_exist():
    tmalign = TMAlign()

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs(['1ycr', '2d7t', '3eiy', '4LWV'])

    chains = {'1ycr': 'A', '2d7t': 'H', '3eiy': 'A', '4LWV': 'A'}

    assert_raises(ValueError, tmalign.tmalign_in_stack, ppdb_stack, chains, target='random')
