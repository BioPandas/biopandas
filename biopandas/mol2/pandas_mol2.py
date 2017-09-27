""" Class for working with Tripos MOL2 files"""
# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import pandas as pd
import numpy as np
from .mol2_io import split_multimol2


COLUMN_NAMES = (
 'atom_id',
 'atom_name',
 'x',
 'y',
 'z',
 'atom_type',
 'subst_id',
 'subst_name',
 'charge'
)

COLUMN_TYPES = (int, str, float, float, float, str, int, str, float)


class PandasMol2(object):
    """
    Object for working with Tripos Mol2 structure files.

    Attributes
    ----------
    df : pandas.DataFrame
        DataFrame of a Mol2's ATOM section

    mol2_text : str
        Mol2 file contents in string format

    code : str
        ID, code, or name of the molecule stored

    pdb_path : str
        Location of the MOL2 file that was read in via `read_mol2`

    """
    def __init__(self):
        self._df = None
        self.mol2_text = ''
        self.header = ''
        self.code = ''
        self.mol2_path = ''

    @property
    def df(self):
        """Acccesses the pandas DataFrame"""
        return self._df

    @df.setter
    def df(self, value):
        """Assign a new value to the pandas DataFrame"""
        raise AttributeError('Please use `PandasMol2._df = ... ` instead\n'
                             'of `PandasMol2.df = ... ` if you are sure that\n'
                             'you want to overwrite the `df` attribute.')
        # self._df = value

    def _load_mol2(self, mol2_lines, mol2_code, columns):
        """Load mol2 contents into assert_raise_message instance"""
        if columns is None:
            col_names = COLUMN_NAMES
            col_types = COLUMN_TYPES
        else:
            col_names, col_types = [], []
            for i in range(len(columns)):
                col_names.append(columns[i][0])
                col_types.append(columns[i][1])

        try:
            self.mol2_text = ''.join(mol2_lines)
            self.code = mol2_code
        except TypeError:
            mol2_lines = [m.decode() for m in mol2_lines]
            self.mol2_text = ''.join(mol2_lines)
            self.code = mol2_code.decode()

        self._df = self._construct_df(mol2_lines, col_names, col_types)

    def read_mol2(self, path, columns=None):
        """Reads Mol2 files (unzipped or gzipped) from local drive

        Note that if your mol2 file contains more than one molecule,
        only the first molecule is loaded into the DataFrame

        Attributes
        ----------
        path : str
            Path to the Mol2 file in .mol2 format or gzipped format (.mol2.gz)

        columns : dict or None (default: None)
            If None, this methods expects a 9-column ATOM section that contains
            the following columns:

            {0:('atom_id', int), 1:('atom_name', str),
             2:('x', float), 3:('y', float), 4:('z', float),
             5:('atom_type', str), 6:('subst_id', int),
             7:('subst_name', str), 8:('charge', float)}

            If your Mol2 files are formatted differently, you can provide your
            own column_mapping dictionary in a format similar to the one above.
            However, note that not all assert_raise_message methods
            may be supported then.

        Returns
        ---------
        self

        """
        mol2_code, mol2_lines = next(split_multimol2(path))
        self._load_mol2(mol2_lines, mol2_code, columns)
        self.mol2_path = path
        return self

    def read_mol2_from_list(self, mol2_lines, mol2_code, columns=None):
        r"""Reads Mol2 file from a list into DataFrames

        Attributes
        ----------
        mol2_lines : list
            A list of lines containing the mol2 file contents. For example,
            ['@<TRIPOS>MOLECULE\n',
             'ZINC38611810\n',
             '   65    68     0     0     0\n',
             'SMALL\n',
             'NO_CHARGES\n',
             '\n',
             '@<TRIPOS>ATOM\n',
             '      1 C1  -1.1786  2.7011  -4.0323 C.3  1 <0>   -0.1537\n',
             '      2 C2  -1.2950  1.2442  -3.5798 C.3  1 <0>   -0.1156\n',
             ...]

        mol2_code : str or None
            Name or ID of the molecule.

        columns : dict or None (default: None)
            If None, this methods expects a 9-column ATOM section that contains
            the following columns:
            {0:('atom_id', int), 1:('atom_name', str),
             2:('x', float), 3:('y', float), 4:('z', float),
             5:('atom_type', str), 6:('subst_id', int),
             7:('subst_name', str), 8:('charge', float)}
            If your Mol2 files are formatted differently, you can provide your
            own column_mapping dictionary in a format similar to the one above.
            However, note that not all assert_raise_message methods may be
            supported then.

        Returns
        ---------
        self

        """
        self._load_mol2(mol2_lines, mol2_code, columns)
        return self

    def _construct_df(self, mol2_lines, col_names, col_types):
        """Construct DataFrames from list of PDB lines."""
        return self._atomsection_to_pandas(self._get_atomsection(mol2_lines),
                                           col_names=col_names,
                                           col_types=col_types)

    @staticmethod
    def _get_atomsection(mol2_lst):
        """Returns atom section from mol2 provided as list of strings"""
        started = False
        for idx, s in enumerate(mol2_lst):
            if s.startswith('@<TRIPOS>ATOM'):
                first_idx = idx + 1
                started = True
            elif started and s.startswith('@<TRIPOS>'):
                last_idx_plus1 = idx
                break
        return mol2_lst[first_idx:last_idx_plus1]

    @staticmethod
    def _atomsection_to_pandas(mol2_atom_lst, col_names, col_types):

        df = pd.DataFrame([lst.split() for lst in mol2_atom_lst],
                          columns=col_names)

        for i in range(df.shape[1]):
            df[col_names[i]] = df[col_names[i]].astype(col_types[i])

        return df

    @staticmethod
    def rmsd(df1, df2, heavy_only=True):
        """Compute the Root Mean Square Deviation between molecules

        Parameters
        ----------
        df1 : pandas.DataFrame
            DataFrame with HETATM, ATOM, and/or ANISOU entries
        df2 : pandas.DataFrame
            Second DataFrame for RMSD computation against df1. Must have the
            same number of entries as df1
        heavy_only : bool (default: True)
            Which atoms to compare to compute the RMSD. If `True` (default),
            computes the RMSD between non-hydrogen atoms only.

        Returns
        ---------
        rmsd : float
            Root Mean Square Deviation between df1 and df2

        """
        if df1.shape[0] != df2.shape[0]:
            raise AttributeError('DataFrames have unequal lengths')

        if heavy_only:
            d1 = df1[df1['atom_type'] != 'H']
            d2 = df2[df2['atom_type'] != 'H']
        else:
            d1, d2 = df1, df2

        total = ((d1['x'].values - d2['x'].values)**2 +
                 (d1['y'].values - d2['y'].values)**2 +
                 (d1['z'].values - d2['z'].values)**2)
        rmsd = round((total.sum() / df1.shape[0])**0.5, 4)
        return rmsd

    def distance(self, xyz=(0.00, 0.00, 0.00)):
        """Computes Euclidean distance between atoms in
            self.df and a 3D point.

        Parameters
        ----------
        xyz : tuple (0.00, 0.00, 0.00)
            X, Y, and Z coordinate of the reference center for the distance
            computation

        Returns
        ---------
        pandas.Series : Pandas Series object containing the Euclidean
            distance between the atoms in the atom section and `xyz`.

        """
        return np.sqrt(np.sum(self.df[['x', 'y', 'z']]
                       .subtract(xyz, axis=1)**2, axis=1))

    @staticmethod
    def distance_df(df, xyz=(0.00, 0.00, 0.00)):
        """Computes Euclidean distance between atoms and a 3D point.

        Parameters
        ----------
        df : DataFrame
            DataFrame containing entries similar to the PandasMol2.df
            format for the
            the distance computation to the `xyz` reference coordinates.
        xyz : tuple (0.00, 0.00, 0.00)
            X, Y, and Z coordinate of the reference center for the distance
            computation

        Returns
        ---------
        pandas.Series : Pandas Series object containing the Euclidean
            distance between the atoms in the atom section and `xyz`.

        """

        return np.sqrt(np.sum(df[['x', 'y', 'z']]
                       .subtract(xyz, axis=1)**2, axis=1))
