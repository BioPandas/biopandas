# Authors: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause

import pandas as pd
import numpy as np
import sys
import gzip
import urllib
from .engines import pdb_records

class PandasPDB(object):
    def __init__(self):
        self._df = {}
        self.pdb_text = ''
        self._get_dict = {}

    @property
    def df(self):
        return self._df

    @property
    def df_lean(self):
        self._df[np.isfinite(self._df)]
        return self._df

    def read_pdb(self, path):
        self.pdb_text = self._read_pdb(path=path)
        self._df = self._construct_df(pdb_lines=self.pdb_text.split('\n'))

    def fetch_pdb(self, pdb_code):
        """Fetches PDB file contents from rcsb.org.

        Parameters
        ----------
        pdb_code : str
            A 4-letter PDB code, e.g., "3eiy"

        """
        self.pdb_text = self._fetch_pdb(pdb_code)
        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines())

    def get(s, df=None):
        if not self._getdict:
            self._init_get_dict()
        if s not in self._getdict.keys():
            raise ValueError('s must be in %s' % self._getdict.keys())
        if not df:
            df = self._df['ATOM']

    def _init_get_dict(self):
        get_dict = {'main-chain': self._get_mainchain,
                    'hydrogens': self._get_hydrogens,
                    'waters': self._get_waters,
                    'c-alpha': self._get_calpha}
        self._getdict = get_dict

    @staticmethod
    def rmsd(df1, df2):
        pass

    @staticmethod
    def _read_pdb(path):
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
        txt = None
        try:
            response = urllib.request.urlopen('http://www.rcsb.org/pdb/files/%s.pdb' % pdb_code.lower())
            txt = response.read().decode('utf-8')
        except urllib.request.HTTPError as e:
            print('HTTP Error %s' %e.code)
        except urllib.request.URLError as e:
            print('URL Error %s' %e.args)
        return txt

    @staticmethod
    def _get_mainchain(df):
        mc =  df[(ppdb.df['ATOM']['atom_name'] == 'C') |
                 (ppdb.df['ATOM']['atom_name'] == 'O') |
                 (ppdb.df['ATOM']['atom_name'] == 'N') |
                 (ppdb.df['ATOM']['atom_name'] == 'CA')]
        return mc

    @staticmethod
    def _get_pdbcode(df):
        pass

    @staticmethod
    def _get_header(df):
        pass

    @staticmethod
    def _get_title(df):
        pass

    @staticmethod
    def _get_moleculename(df):
        pass

    @staticmethod
    def _get_waters(df):
        pass

    @staticmethod
    def _get_hydrogens(df):
        pass

    @staticmethod
    def _get_calpha(df):
        return df[df['atom_name'] == 'CA']

    @staticmethod
    def _construct_df(pdb_lines):
        valids = tuple(pdb_records.keys())
        line_lists = {r:[] for r in valids}
        line_lists['OTHERS'] = []
        for line_num, line in enumerate(pdb_lines):
            if line.strip() and line.startswith(valids):
                record = line[:6].rstrip()
                line_ele = ['' for _ in range(len(pdb_records[record])+1)]
                for idx, ele in enumerate(pdb_records[record]):
                    line_ele[idx] = line[ele['line'][0]:ele['line'][1]].strip()
                line_ele[-1] = line_num
                line_lists[record].append(line_ele)
            else:
                line_lists['OTHERS'].append([line[:6], line[6:-2].strip(), line_num])

        dfs = {}
        for r in line_lists.items():
            df = pd.DataFrame(r[1], columns=[c['id'] for c in pdb_records[r[0]]]+['line_idx'])
            for c in pdb_records[r[0]]:
                try:
                    df[c['id']] = df[c['id']].astype(c['type'])
                except ValueError:
                    # Value Error expected if float/int columns are empty strings
                    df[c['id']] = pd.Series(np.nan, index=df.index)
                    pass
            dfs[r[0]] = df
        return dfs

    def to_pdb(self, path, to_write=None, gz=False, append_newline=False):
        if gz:
            openf = gzip.open
            w_mode = 'wb'
        else:
            openf=open
            w_mode = 'w'
        if not to_write:
            to_write = self._df
            records = self._df.keys()
            dfs = {r:to_write[r].copy() for r in records if not to_write[r].empty}
            for r in dfs.keys():
                for col in pdb_records[r]:
                    dfs[r][col['id']] = dfs[r][col['id']].apply(col['strf'])
                    dfs[r]['OUT'] = pd.Series('', index=dfs[r].index)

                for c in (c for c in dfs[r].columns if c not in {'line_idx', 'OUT'}):
                    dfs[r]['OUT'] = dfs[r]['OUT'] + dfs[r][c]
            df = pd.concat(dfs)
            df.sort(columns='line_idx', inplace=True)
            with openf(path, w_mode) as f:
                f.write('\n'.join(df['OUT'].tolist()))
                if append_newline:
                    f.write('\n')
        """
        if isinstance(to_write, pd.core.frame.DataFrame):
            records = np.unique(to_write['record_name'])
            dfs = [to_write[to_write['record_name'] == r ].copy() for r in records]
        else:
            records = to_write
            dfs = [self._df[r].copy() for r in records]
        with open(path, 'w') as f:
            for df, r in zip(dfs, records):
                col_to_write = []
                for col in pdb_records[r]:
                    df[col['id']] = df[col['id']].apply(col['strf'])
                    col_to_write.append(col['id'])
                for ele in df[col_to_write].values:
                    f.write(''.join(ele) + '\n')
        """
