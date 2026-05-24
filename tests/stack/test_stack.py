# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
 # Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

from biopandas.stack.stack import PandasPdbStack
from biopandas.align import TMAlign
import os
import pytest

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), "data", "3eiy.pdb")

TESTDATA_FILENAME3 = os.path.join(
    os.path.dirname(__file__), "data", "1ycr.pdb"
)
TESTDATA_FILENAME4 = os.path.join(
    os.path.dirname(__file__), "data", "2d7t.pdb"
)

def has_usalign():
    # Check if USalign exists in expected locations
    paths = [
        './biopandas/align/USalign',
        './biopandas/align/USalign.exe',
        '../../biopandas/align/USalign',
        '../../biopandas/align/USalign.exe'
    ]
    return any(os.path.exists(p) for p in paths)

# Helper functions for testing
def filter_by_chain(key, pdb, chain):
    # Example function for applying filtering
    pdb.df['ATOM'] = pdb.df['ATOM'].query('chain_id==@chain')
    return pdb

def calculate_chain_lengths(key, pdb):
    # Assuming `pdb` is a PandasPdb object
    lengths = {}
    for ch in pdb.df['ATOM']['chain_id'].unique():
      ch_len = len(pdb.df['ATOM'].query('chain_id==@ch and atom_name=="CA"'))
      lengths[ch] = ch_len
    return lengths


def test_apply_filter():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, '1YCR', 'P99999'])

    args = {'chain': 'A'}
    filtered_stack = stack.apply_filter(filter_by_chain, keep_null=True, **args)
    assert len(filtered_stack.pdbs) == 3

def test_apply_filter_no_keep_null():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, '1YCR', 'P99999'])

    args = {'chain': 'B'}

    filtered_stack = stack.apply_filter(filter_by_chain, keep_null=False, **args)
    assert len(filtered_stack.pdbs) == 1

def test_apply_calculation():
    stack = PandasPdbStack()
    stack.add_pdbs(['1YCR', '1A2B'])

    chain_lengths = stack.apply_calculation(calculate_chain_lengths)
    assert chain_lengths == {'1YCR': {'A': 85, 'B': 13}, '1A2B': {'A': 178}}

def test_delete_entry():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, '1YCR', 'P99999'])
    stack.delete_entry('1YCR')
    assert len(stack.pdbs) == 2
    assert '1YCR' not in stack.pdbs

def test_delete_entry_nonexistent():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, '1YCR', 'P99999'])
    stack.delete_entry('1YCR')
    assert len(stack.pdbs) == 2
    assert '1YCR' not in stack.pdbs
    stack.delete_entry('1YCR')
    assert len(stack.pdbs) == 2
    assert '1YCR' not in stack.pdbs

def test_update_entry():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, '1YCR', 'P99999'])
    stack.update_entry('1YCR', '1A2B')
    assert len(stack.pdbs) == 3
    assert '1YCR' in stack.pdbs
    assert '1A2B' not in stack.pdbs

def test_update_entry_nonexistent():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, '1YCR', 'P99999'])
    stack.update_entry('1A2B', '1YCR')
    assert len(stack.pdbs) == 4
    assert '1YCR' in stack.pdbs
    assert '1A2B' in stack.pdbs


def filter_by_chains(key, pdb, chains):
    # Example function for applying filtering
    chain = chains[key]
    pdb.df['ATOM'] = pdb.df['ATOM'].query('chain_id==@chain')
    return pdb

@pytest.mark.skipif(not has_usalign(), reason="USalign executable not found")
def test_tmalign_inside_multiple_chains():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME3, TESTDATA_FILENAME4])
    stack.add_pdb(TESTDATA_FILENAME3, '1ycr_copy')

    args = {'chains': {'1ycr': 'A', '2d7t': 'H', '3eiy': 'A', '1ycr_copy': 'A'}}
    filtered_stack = stack.apply_filter(filter_by_chains, keep_null=False, **args)
    chains_lens_filtered = filtered_stack.apply_calculation(calculate_chain_lengths)
    assert len(filtered_stack.pdbs) == 4

    tmalign = TMAlign()
    _, transformed_structures, tm_scores = tmalign.tmalign_in_stack(filtered_stack, mobile_chains=args['chains'])

    assert tm_scores['3eiy'] == 0.23483
    assert tm_scores['2d7t'] == 0.27812
    assert tm_scores['1ycr_copy'] == 1

@pytest.mark.skipif(not has_usalign(), reason="USalign executable not found")
def test_tmalign_inside_multiple_chains_specific_target():
    stack = PandasPdbStack()
    stack.add_pdbs([TESTDATA_FILENAME, TESTDATA_FILENAME3, TESTDATA_FILENAME4])
    stack.add_pdb(TESTDATA_FILENAME, '3eiy_copy')

    args = {'chains': {'1ycr': 'A', '2d7t': 'H', '3eiy': 'A', '3eiy_copy': 'A'}}
    filtered_stack = stack.apply_filter(filter_by_chains, keep_null=False, **args)
    chains_lens_filtered = filtered_stack.apply_calculation(calculate_chain_lengths)
    assert len(filtered_stack.pdbs) == 4

    tmalign = TMAlign()
    _, transformed_structures, tm_scores = tmalign.tmalign_in_stack(filtered_stack, target='3eiy', mobile_chains=args['chains'])

    assert tm_scores['1ycr'] == 0.37341
    assert tm_scores['2d7t'] == 0.32871
    assert tm_scores['3eiy_copy'] == 1