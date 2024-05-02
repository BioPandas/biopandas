# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.pdb import PandasPdb
from biopandas.stack.stack import PandasPdbStack
from biopandas.align import TMAlign
import os
import sys
from nose.tools import assert_raises

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
        tmalign = TMAlign(os.path.join(os.path.dirname(__file__), "../USalign.exe"))
    else:
        tmalign = TMAlign(os.path.join(os.path.dirname(__file__), "../USalign"))

    assert 'USalign' in tmalign.tmalign_path

def test_init_with_path_not_ok():
    assert_raises(FileNotFoundError, TMAlign, 'fake_Usalign')

def test_init_with_path_none():
    os.chdir('..')

    if os.path.exists(os.path.join(os.path.dirname(__file__), "../USalign.exe")):
        tmalign = TMAlign()
        assert 'USalign' in tmalign.tmalign_path
    elif os.path.exists(os.path.join(os.path.dirname(__file__), "../USalign")):
        tmalign = TMAlign()
        assert 'USalign' in tmalign.tmalign_path
    else:
        assert_raises(ValueError, TMAlign, None)

def test_parse_tmalign_rotation_matrix():
    tmalign = TMAlign()
    matrix, translation = tmalign.parse_tmalign_rotation_matrix(
        os.path.join(os.path.dirname(__file__), "data", "tmalign_output.txt")
    )
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

    target_file = tmalign.write_pdb_to_temp_file(ppdb)
    transformed_mobile, tm_score = tmalign.process_structure_for_tmalign(target_file, ppdb2)

    assert tm_score == 1

def test_process_structure_for_tmalign_not_perfect():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)
    ppdb2 = PandasPdb()
    ppdb2.read_pdb(TESTDATA_FILENAME3)

    target_file = tmalign.write_pdb_to_temp_file(ppdb)
    transformed_mobile, tm_score = tmalign.process_structure_for_tmalign(target_file, ppdb2)

    assert tm_score == 0.07349

def test_tmalign_to_stack_chains_exist():
    tmalign = TMAlign()
    ppdb = PandasPdb()
    ppdb.read_pdb(TESTDATA_FILENAME)

    ppdb_stack = PandasPdbStack()
    ppdb_stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME3, TESTDATA_FILENAME4])

    transformed_structures, tm_scores = tmalign.tmalign_to(ppdb, ppdb_stack,
                                                           'A',
                                                           {'2jyf': 'A', '2d7t': 'H', '3eiy': 'A'})
    assert tm_scores['2jyf'] == 0.07349
    assert tm_scores['2d7t'] == 0.24401
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
    assert tm_score == 0.07349

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