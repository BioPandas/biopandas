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
        """Write a PandasPdb object's data to a temporary PDB file and return the file handle.
        :param pdb: the PandasPdb object to write to the file.

        :return: file handle
        """
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        pdb.to_pdb(path=temp_file.name, records=None, gz=False, append_newline=True)
        return temp_file

    def filter_and_validate_chain(self, pdb, chain_id):
        """Filter the PandasPdb by chain_id and validate the presence of the chain.
        :param pdb: the PandasPdb object to filter.
        :param chain_id: the chain ID to filter by.

        :return: filtered_pdb
        """
        filtered_pdb = deepcopy(pdb)
        filtered_atoms = pdb.df['ATOM'][pdb.df['ATOM']['chain_id'].isin([chain_id])]
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
