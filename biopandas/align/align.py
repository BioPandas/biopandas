""" Class for aligning PDB structures"""

# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import tempfile
from copy import deepcopy

import numpy as np


class Align():
    def __init__(self):
        pass

    def write_pdb_to_temp_file(self, pdb):
        """Write a PandasPdb/PandasMmcif object's data to a temporary structure file and return the file handle.
        :param pdb: the Pandas object to write to the file.

        :return: file handle
        """
        # if pdb is PandasPdb object,call to_pdb, if PandasMmcif, call to_mmcif
        if hasattr(pdb, 'to_pdb'):
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
            pdb.to_pdb(path=temp_file.name, records=None, gz=False, append_newline=True)
        elif hasattr(pdb, 'to_mmcif'):
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.cif')
            pdb.to_mmcif(path=temp_file.name, records=None, gz=False)

        return temp_file

    def filter_and_validate_chain(self, pdb, chain_id):
        """Filter the PandasPdb by chain_id and validate the presence of the chain.
        :param pdb: the PandasPdb object to filter.
        :param chain_id: the chain ID to filter by.

        :return: filtered_pdb
        """
        # need to find if chain is called chain_id, label_asym_id or auth_asym_id.
        # order for checking: chain_id, auth_asym_id, label_asym_id
        chain_col = None
        if 'chain_id' not in pdb.df['ATOM'].columns:
            if 'auth_asym_id' in pdb.df['ATOM'].columns:
                chain_col = 'auth_asym_id'
            elif 'label_asym_id' in pdb.df['ATOM'].columns:
                chain_col = 'label_asym_id'
            else:
                raise ValueError("No recognized chain identifier column found in the ATOM dataframe.")
        else:
            chain_col = 'chain_id'

        filtered_pdb = deepcopy(pdb)
        filtered_atoms = pdb.df['ATOM'][pdb.df['ATOM'][chain_col].isin([chain_id])]
        if filtered_atoms.empty:
            raise ValueError(f"No such chain '{chain_id}' found in the structure.")
        filtered_pdb.df['ATOM'] = filtered_atoms
        return filtered_pdb

    def transform(self, coords, matrix, translation):
        """Apply the rotation matrix and translation vector to the structure.
        :param coords: the coordinates to transform.
        :param matrix: the rotation matrix.
        :param translation: the translation vector.

        :return: transformed coordinates as a numpy array.
        """

        return np.dot(coords, matrix.T) + translation
