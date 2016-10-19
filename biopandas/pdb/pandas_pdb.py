# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import pandas as pd
import numpy as np
import sys
import gzip
from warnings import warn
try:
    from urllib.request import urlopen
    from urllib.error import HTTPError, URLError
except ImportError:
    from urllib2 import urlopen, HTTPError, URLError  # Python 2.7 compatible
from .engines import pdb_records


class PandasPDB(object):
    """ Object for working with Protein Databank structure files.

    Attributes
    ----------
    df : dict
        Dictionary storing pandas DataFrames for PDB record sections.
        The dictionary keys are {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}
        where 'OTHERS' contains all entries that are not parsed as
        'ATOM', 'HETATM', or 'ANISOU'

    pdb_text : str
        PDB file contents in raw text format

    header : str
        PDB file description

    code : str
        PDB code

    """
    def __init__(self):
        self._df = {}
        self.pdb_text = ''
        self.header = ''
        self.code = ''
        self._get_dict = {}

    @property
    def df(self):
        """Acccess dictionary of pandas DataFrames for PDB record sections."""
        return self._df

    def read_pdb(self, path):
        """Read PDB files (unzipped or gzipped) from local drive

        Attributes
        ----------
        path : str
            Path to the PDB file in .pdb format or gzipped format (.pdb.gz)

        Returns
        ---------
        self

        """
        self.pdb_text = self._read_pdb(path=path)
        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines(True))
        self.header, self.code = self._parse_header_code()
        return self

    def fetch_pdb(self, pdb_code):
        """Fetches PDB file contents from the Protein Databank at rcsb.org.

        Parameters
        ----------
        pdb_code : str
            A 4-letter PDB code, e.g., "3eiy"

        Returns
        ---------
        self

        """
        self.pdb_text = self._fetch_pdb(pdb_code)
        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines(True))
        return self

    def get(self, s, df=None, invert=False):
        """Filter PDB DataFrames by properties

        Parameters
        ----------
        s : str  in {'main chain', 'hydrogen', 'c-alpha', 'heavy'}
            String to specify which entries to return

        df : pandas.DataFrame, default: None
            Optional DataFrame to perform the filter operation on.
            If df=None, filters on self.df['ATOM']

        invert : bool, default: True
            Inverts the search query. For example if s='hydrogen' and
            invert=True, all but hydrogen entries are returned

        Returns
        --------
        df : pandas.DataFrame
            Returns a DataFrame view on the filtered entries.

        """
        if not self._get_dict:
            self._get_dict = self._init_get_dict()
        if s not in self._get_dict.keys():
            raise AttributeError('s must be in %s' % self._get_dict.keys())
        if not df:
            df = self._df['ATOM']
        return self._get_dict[s](df, invert=invert)

    def impute_element(self, sections=['ATOM', 'HETATM'], inplace=False):
        """Impute element_symbol from atom_name section.

        Parameters
        ----------
        sections : iterable (default: ['ATOM', 'HETATM'])
            Coordinate sections for which the element symbols should be
            imputed.

        inplace : bool (default: False)
            Performs the operation in-place if True and returns a copy of the
            PDB DataFrame otherwise.

        Returns
        ---------
        DataFrame

        """
        if inplace:
            t = self.df
        else:
            t = self.df.copy()
            for d in self.df:
                t[d] = self.df[d].copy()

        for sec in sections:
            t[sec]['element_symbol'] = \
                t[sec][['atom_name', 'element_symbol']].\
                apply(lambda x: x[0][1]
                      if len(x[1]) == 3
                      else x[0][0], axis=1)
        return t

    @staticmethod
    def rmsd(df1, df2, s=None, invert=False):
        """Compute the Root Mean Square Deviation between molecules.

        Parameters
        ----------
        df1 : pandas.DataFrame
            DataFrame with HETATM, ATOM, and/or ANISOU entries

        df2 : pandas.DataFrame
            Second DataFrame for RMSD computation against df1. Must have the
            same number of entries as df1

        s : {'main chain', 'hydrogen', 'c-alpha', 'heavy', 'carbon'} or None,
            default: None
            String to specify which entries to consider. If None, considers
            all atoms for comparison.

        invert : bool, default: False
            Inverts the string query if true. For example, the setting
            `s='hydrogen', invert=True` computes the RMSD based on all
            but hydrogen atoms.

        Returns
        ---------
        rmsd : float
            Root Mean Square Deviation between df1 and df2

        """
        if df1.shape[0] != df2.shape[0]:
            raise AttributeError('DataFrames have unequal lengths')
        get_dict = PandasPDB._init_get_dict()
        if s:
            if s not in get_dict.keys():
                raise AttributeError('s must be in '
                                     '%s or None' % get_dict.keys())
            df1 = get_dict[s](df1, invert=invert)
            df2 = get_dict[s](df2, invert=invert)

        total = ((df1['x_coord'] - df2['x_coord'])**2 +
                 (df1['y_coord'] - df2['y_coord'])**2 +
                 (df1['z_coord'] - df2['z_coord'])**2)
        rmsd = round((total.sum() / df1.shape[0])**0.5, 4)
        return rmsd

    @staticmethod
    def _init_get_dict():
        """Initialize dictionary for filter operations."""
        get_dict = {'main chain': PandasPDB._get_mainchain,
                    'hydrogen': PandasPDB._get_hydrogen,
                    'c-alpha': PandasPDB._get_calpha,
                    'carbon': PandasPDB._get_carbon,
                    'heavy': PandasPDB._get_heavy}
        return get_dict

    @staticmethod
    def _read_pdb(path):
        """Read PDB file from local drive."""
        r_mode = 'r'
        openf = open
        if path.endswith('.gz'):
            r_mode = 'rb'
            openf = gzip.open
        with openf(path, r_mode) as f:
            txt = f.read()
        if path.endswith('.gz'):
            if sys.version_info[0] >= 3:
                txt = txt.decode('utf-8')
            else:
                txt = txt.encode('ascii')
        return txt

    @staticmethod
    def _fetch_pdb(pdb_code):
        """Load PDB file from rcsb.org."""
        txt = None
        try:
            response = urlopen('http://www.rcsb.org/pdb/files/%s.pdb'
                               % pdb_code.lower())
            txt = response.read()
            if sys.version_info[0] >= 3:
                txt = txt.decode('utf-8')
            else:
                txt = txt.encode('ascii')
        except HTTPError as e:
            print('HTTP Error %s' % e.code)
        except URLError as e:
            print('URL Error %s' % e.args)
        return txt

    def _parse_header_code(self):
        """Extract header information and PDB code."""
        code, header = '', ''
        if 'OTHERS' in self.df:

            header = (self.df['OTHERS'][self.df['OTHERS']['record_name'] ==
                      'HEADER'])
            if not header.empty:
                header = header['entry'].values[0]
                s = header.split()
                if s:
                    code = s[-1].lower()
        return header, code

    @staticmethod
    def _get_mainchain(df, invert):
        """Return only main chain atom entries from a DataFrame"""
        if invert:
            mc = df[(df['atom_name'] != 'C') &
                    (df['atom_name'] != 'O') &
                    (df['atom_name'] != 'N') &
                    (df['atom_name'] != 'CA')]
        else:
            mc = df[(df['atom_name'] == 'C') |
                    (df['atom_name'] == 'O') |
                    (df['atom_name'] == 'N') |
                    (df['atom_name'] == 'CA')]
        return mc

    @staticmethod
    def _get_hydrogen(df, invert):
        """Return only hydrogen atom entries from a DataFrame"""
        if invert:
            return df[(df['element_symbol'] != 'H')]
        else:
            return df[(df['element_symbol'] == 'H')]

    @staticmethod
    def _get_heavy(df, invert):
        """Return only heavy atom entries from a DataFrame"""
        if invert:
            return df[df['element_symbol'] == 'H']
        else:
            return df[df['element_symbol'] != 'H']

    @staticmethod
    def _get_calpha(df, invert):
        """Return c-alpha atom entries from a DataFrame"""
        if invert:
            return df[df['atom_name'] != 'CA']
        else:
            return df[df['atom_name'] == 'CA']

    @staticmethod
    def _get_carbon(df, invert):
        """Return c-alpha atom entries from a DataFrame"""
        if invert:
            return df[df['element_symbol'] == 'C']
        else:
            return df[df['element_symbol'] != 'C']

    @staticmethod
    def _construct_df(pdb_lines):
        """Construct DataFrames from list of PDB lines."""
        valids = tuple(pdb_records.keys())
        line_lists = {r: [] for r in valids}
        line_lists['OTHERS'] = []
        for line_num, line in enumerate(pdb_lines):
            if line.strip():
                if line.startswith(valids):
                    record = line[:6].rstrip()
                    line_ele = ['' for _ in range(len(
                        pdb_records[record]) + 1)]
                    for idx, ele in enumerate(pdb_records[record]):
                        line_ele[idx] = (line[ele['line'][0]:ele['line'][1]]
                                         .strip())
                    line_ele[-1] = line_num
                    line_lists[record].append(line_ele)
                else:
                    line_lists['OTHERS'].append([line[:6].rstrip(),
                                                line[6:-1].rstrip(), line_num])

        dfs = {}
        for r in line_lists.items():
            df = pd.DataFrame(r[1], columns=[c['id'] for c in
                                             pdb_records[r[0]]] + ['line_idx'])
            for c in pdb_records[r[0]]:
                try:
                    df[c['id']] = df[c['id']].astype(c['type'])
                except ValueError:
                    # expect ValueError if float/int columns are empty strings
                    df[c['id']] = pd.Series(np.nan, index=df.index)

            dfs[r[0]] = df
        return dfs

    def to_pdb(self, path, records=None, gz=False, append_newline=True):
        """Write record DataFrames to a PDB file or gzipped PDB file.

        Parameters
        ----------
        path : str
            A valid output path for the pdb file

        records : iterable, default: None
            A list of PDB record sections in
            {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'} that are to be written.
            Writes all lines to PDB if records=None

        gz : bool, default: False
            Writes a gzipped PDB file if True

        append_newline : bool, default: True
            Appends a new line at the end of the PDB file if True

        """
        if gz:
            openf = gzip.open
            w_mode = 'wt'
        else:
            openf = open
            w_mode = 'w'
        if not records:
            records = self.df.keys()

        dfs = {r: self.df[r].copy() for r in records if not self.df[r].empty}

        for r in dfs.keys():
            for col in pdb_records[r]:
                dfs[r][col['id']] = dfs[r][col['id']].apply(col['strf'])
                dfs[r]['OUT'] = pd.Series('', index=dfs[r].index)

            for c in (c for c in dfs[r].columns
                      if c not in {'line_idx', 'OUT'}):
                dfs[r]['OUT'] = dfs[r]['OUT'] + dfs[r][c]

        df = pd.concat(dfs)
        if pd.__version__ < '0.17.0':
            warn("You are using an old pandas version (< 0.17)"
                 " that relies on the old sorting syntax."
                 " Please consider updating your pandas"
                 " installation to a more recent version.",
                 DeprecationWarning)
            df.sort(columns='line_idx', inplace=True)
        else:
            df.sort_values(by='line_idx', inplace=True)
        with openf(path, w_mode) as f:

            s = df['OUT'].tolist()
            for idx in range(len(s)):
                if len(s[idx]) < 80:
                    s[idx] = '%s%s' % (s[idx], ' ' * (80 - len(s[idx])))
            to_write = '\n'.join(s)
            f.write(to_write)
            if append_newline:
                if gz:
                    f.write('\n')
                else:
                    f.write('\n')
