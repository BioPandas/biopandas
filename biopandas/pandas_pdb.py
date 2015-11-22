# Authors: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause

import pandas as pd
import numpy as np
import sys
import gzip
try:
    from urllib.request import urlopen
    from urllib.error import HTTPError, URLError
except ImportError:
    from urllib2 import urlopen, Request, HTTPError, URLError
from .engines import pdb_records

class PandasPDB(object):
    def __init__(self):
        self._df = {}
        self.pdb_text = ''
        self.title = ''
        self.code = ''
        self._get_dict = {}

    @property
    def df(self):
        return self._df

    def read_pdb(self, path):
        self.pdb_text = self._read_pdb(path=path)
        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines(True))
        self.title, self.code = self._parse_title_code()

    def fetch_pdb(self, pdb_code):
        """Fetches PDB file contents from rcsb.org.

        Parameters
        ----------
        pdb_code : str
            A 4-letter PDB code, e.g., "3eiy"

        """
        self.pdb_text = self._fetch_pdb(pdb_code)
        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines(True))

    def get(self, s, df=None):
        if not self._get_dict:
            self._get_dict = self._init_get_dict()
        if s not in self._get_dict.keys():
            raise AttributeError('s must be in %s' % self._get_dict.keys())
        if not df:
            df = self._df['ATOM']
        return self._get_dict[s](df)

    @staticmethod
    def rmsd(df1, df2, s='no hydrogen'):
        if  df1.shape[0] != df2.shape[0]:
            raise AttributeError('DataFrames have unequal lengths')
        get_dict = PandasPDB._init_get_dict()
        if s:
            if s not in get_dict.keys():
                raise AttributeError('s must be in %s or None' % get_dict.keys())
            df1 = get_dict[s](df1)
            df2 = get_dict[s](df2)

        total = ((df1['x_coord'] - df2['x_coord'])**2
               + (df1['y_coord'] - df2['y_coord'])**2
               + (df1['z_coord'] - df2['z_coord'])**2)
        rmsd = round(( total.sum() / df1.shape[0] )**0.5, 4)
        return rmsd


    @staticmethod
    def _init_get_dict():
        get_dict = {'main chain': PandasPDB._get_mainchain,
                    'hydrogen': PandasPDB._get_hydrogen,
                    'no hydrogen': PandasPDB._get_no_hydrogen,
                    'c-alpha': PandasPDB._get_calpha}
        return get_dict

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
            response = urlopen('http://www.rcsb.org/pdb/files/%s.pdb' % pdb_code.lower())
            txt = response.read()
            if sys.version_info[0] >= 3:
                txt = txt.decode('utf-8')
            else:
                txt = txt.encode('ascii')
        except HTTPError as e:
            print('HTTP Error %s' %e.code)
        except URLError as e:
            print('URL Error %s' %e.args)
        return txt

    def _parse_title_code(self):
        code, title = '', ''
        if 'OTHERS' in self.df:

            header = self.df['OTHERS'][self.df['OTHERS']['record_name'] == 'HEADER']
            if not header.empty:
                title = header['entry'].values[0]
                s = title.split()
                if s:
                    code = s[-1].lower()
        return title, code


    @staticmethod
    def _get_mainchain(df):
        mc =  df[(df['atom_name'] == 'C') |
                 (df['atom_name'] == 'O') |
                 (df['atom_name'] == 'N') |
                 (df['atom_name'] == 'CA')]
        return mc


    @staticmethod
    def _get_hydrogen(df):
        df_h = df[(df['atom_name'] == 'H')]
        return df_h

    @staticmethod
    def _get_no_hydrogen(df):
        df_noh = df[(df['atom_name'] != 'H')]
        return df_noh

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
                line_lists['OTHERS'].append([line[:6], line[6:-1].strip(), line_num])

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
