""" Class for working with PDB files"""

# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas
from __future__ import annotations

import gzip
import sys
import textwrap
import warnings
from copy import deepcopy
from io import StringIO
from typing import List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen
from warnings import warn

import numpy as np
import pandas as pd
from looseversion import LooseVersion

from biopandas.constants import ATOMIC_MASSES

from .engines import amino3to1dict, pdb_df_columns, pdb_records

pd_version = LooseVersion(pd.__version__)


class PandasPdb(object):
    """
    Object for working with Protein Databank structure files.

    Attributes
    ----------
    df : dict
        Dictionary storing pandas DataFrames for PDB record sections.
        The dictionary keys are {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'}
        where 'OTHERS' contains all entries that are not parsed as
        'ATOM', 'HETATM', or 'ANISOU'.

    pdb_text : str
        PDB file contents in raw text format.

    pdb_path : Union[str, os.PathLike]
        Location of the PDB file that was read in via `read_pdb`
        or URL of the page where the PDB content was fetched from
        if `fetch_pdb` was called.

    header : str
        PDB file description.

    code : str
        PDB code

    """

    def __init__(self):
        self._df = {}
        self.pdb_text = ""
        self.header = ""
        self.code = ""
        self._get_dict = {}
        self.pdb_path = ""

    @property
    def df(self):
        """Access dictionary of pandas DataFrames for PDB record sections."""
        return self._df

    @df.setter
    def df(self, value):
        """Assign a new value to the pandas DataFrame"""
        raise AttributeError(
            "Please use `PandasPdb._df = ... ` instead\n"
            "of `PandasPdb.df = ... ` if you are sure that\n"
            "you want to overwrite the `df` attribute."
        )
        # self._df = value

    def read_pdb(self, path):
        """Read PDB files (unzipped or gzipped) from local drive

        Attributes
        ----------
        path : str
            Path to the PDB file in .pdb format or gzipped format (.pdb.gz).

        Returns
        ---------
        self

        """
        self.pdb_path, self.pdb_text = self._read_pdb(path=str(path))
        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines(True))
        self.header, self.code = self._parse_header_code()
        return self

    def read_pdb_from_list(self, pdb_lines):
        """Reads PDB file from a list into DataFrames

        Attributes
        ----------
        pdb_lines : list
            A list of lines containing the pdb file contents.

        Returns
        ---------
        self

        """
        self.pdb_text = "".join(pdb_lines)
        self._df = self._construct_df(pdb_lines)
        self.header, self.code = self._parse_header_code()
        return self

    def fetch_pdb(
        self,
        pdb_code: Optional[str] = None,
        uniprot_id: Optional[str] = None,
        source: str = "pdb",
    ):
        """Fetches PDB file contents from the Protein Databank at rcsb.org or AlphaFold database
                at https://alphafold.ebi.ac.uk/.
        .

                Parameters
                ----------
                pdb_code : str, optional
                    A 4-letter PDB code, e.g., `"3eiy"` to retrieve structures from the PDB.
                    Defaults to `None`.

                uniprot_id : str, optional
                    A UniProt Identifier, e.g., `"Q5VSL9"` to retrieve structures from the AF2 database.
                    Defaults to `None`.

                source : str
                    The source to retrieve the structure from
                    (`"pdb"`, `"alphafold2-v3"`, `"alphafold2-v4"`(latest)).
                    Defaults to `"pdb"`.

                Returns
                ---------
                self

        """
        # Sanitize input
        invalid_input_identifier_1 = pdb_code is None and uniprot_id is None
        invalid_input_identifier_2 = pdb_code is not None and uniprot_id is not None
        invalid_input_combination_1 = uniprot_id is not None and source == "pdb"
        invalid_input_combination_2 = pdb_code is not None and source in {
            "alphafold2-v3",
            "alphafold2-v4",
        }

        if invalid_input_identifier_1 or invalid_input_identifier_2:
            raise ValueError("Please provide either a PDB code or a UniProt ID.")

        if invalid_input_combination_1:
            raise ValueError(
                "Please use a 'pdb_code' instead of 'uniprot_id' for source='pdb'."
            )
        elif invalid_input_combination_2:
            raise ValueError(
                f"Please use a 'uniprot_id' instead of 'pdb_code' for source={source}."
            )

        if source == "alphafold2-v3":
            af2_version = 3
            self.pdb_path, self.pdb_text = self._fetch_af2(uniprot_id, af2_version)
        elif source == "alphafold2-v4":
            af2_version = 4
            self.pdb_path, self.pdb_text = self._fetch_af2(uniprot_id, af2_version)
        elif source == "pdb":
            self.pdb_path, self.pdb_text = self._fetch_pdb(pdb_code)
        else:
            raise ValueError(
                f"Invalid source: {source}."
                " Please use one of 'pdb' or 'alphafold2-v3' or 'alphafold2-v4'."
            )

        self._df = self._construct_df(pdb_lines=self.pdb_text.splitlines(True))
        return self

    def get(self, s, df=None, invert=False, records=("ATOM", "HETATM")):
        """Filter PDB DataFrames by properties

        Parameters
        ----------
        s : str  in {'main chain', 'hydrogen', 'c-alpha', 'heavy'}
            String to specify which entries to return.

        df : pandas.DataFrame, default: None
            Optional DataFrame to perform the filter operation on.
            If df=None, filters on self.df['ATOM'].

        invert : bool, default: True
            Inverts the search query. For example if s='hydrogen' and
            invert=True, all but hydrogen entries are returned.

        records : iterable, default: ('ATOM', 'HETATM')
            Specify which record sections to consider. For example, to consider
            both protein and ligand atoms, set `records=('ATOM', 'HETATM')`.
            This setting is ignored if `df` is not set to None.
            For downward compatibility, a string argument is still supported
            but deprecated and will be removed in future versions.

        Returns
        --------
        df : pandas.DataFrame
            Returns a DataFrame view on the filtered entries.

        """
        if isinstance(records, str):
            warnings.warn(
                "Using a string as `records` argument is "
                "deprecated and will not be supported in future"
                " versions. Please use a tuple or"
                " other iterable instead",
                DeprecationWarning,
            )
            records = (records,)

        if not self._get_dict:
            self._get_dict = self._init_get_dict()
        if s not in self._get_dict.keys():
            raise AttributeError(f"s must be in {self._get_dict.keys()}")
        if not df:
            df = pd.concat(objs=[self.df[i] for i in records])
        return self._get_dict[s](df, invert=invert)

    def impute_element(self, records=("ATOM", "HETATM"), inplace=False):
        """Impute element_symbol from atom_name section.

        Parameters
        ----------
        records : iterable, default: ('ATOM', 'HETATM')
            Coordinate sections for which the element symbols should be
            imputed.

        inplace : bool, default: False
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

        for sec in records:
            t[sec]["element_symbol"] = t[sec][["atom_name", "element_symbol"]].apply(
                lambda x: x[0][1] if len(x[1]) == 3 else x[0][0], axis=1
            )
        return t

    def add_remark(self, code, text="", indent=0):
        """Add custom REMARK entry.

        The remark will be inserted to preserve the ordering of REMARK codes, i.e. if the code is
        `n` it will be added after all remarks with codes less or equal to `n`. If the object does
        not store any remarks the remark will be inserted right before the first of ATOM, HETATM or
        ANISOU records.

        Parameters
        ----------
        code : int
            REMARK code according to PDB standards.

        text : str
            The text of the remark. If the text does not fit into a single line it will be wrapped
            into multiple lines of REMARK entries. Likewise, if the text contains new line
            characters it will be split accordingly.

        indent : int, default: 0
            Number of white spaces inserted before the text of the remark.

        Returns
        ---------
        Nothing

        """
        # Prepare info from self
        if "OTHERS" in self.df:
            df_others = self.df["OTHERS"]
        else:
            df_others = pd.DataFrame(columns=["record_name", "entry", "line_idx"])
        record_types = list(
            filter(lambda x: x in self.df, ["ATOM", "HETATM", "ANISOU"])
        )
        remarks = df_others[df_others["record_name"] == "REMARK"]["entry"]

        # Find index and line_idx where to insert the remark to preserve remark code order
        if len(remarks):
            remark_codes = remarks.apply(lambda x: x.split(maxsplit=1)[0]).astype(int)
            insertion_pos = remark_codes.searchsorted(code, side="right")
            if insertion_pos < len(remark_codes):  # Remark in the middle
                insertion_idx = remark_codes.index[insertion_pos]
                insertion_line_idx = df_others.loc[insertion_idx]["line_idx"]
            else:  # Last remark
                insertion_idx = len(remark_codes)
                insertion_line_idx = df_others["line_idx"].iloc[-1] + 1
        else:  # First remark
            insertion_idx = 0
            insertion_line_idx = min(
                [self.df[r]["line_idx"].min() for r in record_types]
            )

        # Wrap remark to fit into 80 characters per line and add indentation
        wrapper = textwrap.TextWrapper(width=80 - (11 + indent))
        lines = sum(
            [wrapper.wrap(line.strip()) or [" "] for line in text.split("\n")], []
        )
        lines = list(map(lambda x: f"{code:4} " + indent * " " + x, lines))

        # Shift data frame indices and row indices to create space for the remark
        # Create space in OTHERS
        line_idx = df_others["line_idx"].copy()
        line_idx[line_idx >= insertion_line_idx] += len(lines)
        df_others["line_idx"] = line_idx
        index = pd.Series(df_others.index.copy())
        index[index >= insertion_idx] += len(lines)
        df_others.index = index
        # Shift all other record types that follow inserted remark
        for records in record_types:
            df_records = self.df[records]
            if not insertion_line_idx > df_records["line_idx"].max():
                df_records["line_idx"] += len(lines)

        # Put remark into 'OTHERS' data frame
        df_remark = {
            idx: ["REMARK", line, line_idx]
            for idx, line, line_idx in zip(
                range(insertion_idx, insertion_idx + len(lines)),
                lines,
                range(insertion_line_idx, insertion_line_idx + len(lines)),
            )
        }
        df_remark = pd.DataFrame.from_dict(
            df_remark, orient="index", columns=df_others.columns
        )
        self.df["OTHERS"] = pd.concat([df_others, df_remark]).sort_index()

    @staticmethod
    def rmsd(df1, df2, s=None, invert=False, decimals=4):
        """Compute the Root Mean Square Deviation between molecules.

        Parameters
        ----------
        df1 : pandas.DataFrame
            DataFrame with HETATM, ATOM, and/or ANISOU entries.

        df2 : pandas.DataFrame
            Second DataFrame for RMSD computation against df1. Must have the
            same number of entries as df1.

        s : {'main chain', 'hydrogen', 'c-alpha', 'heavy', 'carbon'} or None,
            default: None
            String to specify which entries to consider. If None, considers
            all atoms for comparison.

        invert : bool, default: False
            Inverts the string query if true. For example, the setting
            `s='hydrogen', invert=True` computes the RMSD based on all
            but hydrogen atoms.

        decimals : int, default: 4
            Specifies the number of decimal places to round the final value to.

        Returns
        ---------
        rmsd : float
            Root Mean Square Deviation between df1 and df2

        """
        if df1.shape[0] != df2.shape[0]:
            raise AttributeError("DataFrames have unequal lengths")
        get_dict = PandasPdb._init_get_dict()
        if s:
            if s not in get_dict.keys():
                raise AttributeError(f"s must be in {get_dict.keys()} or None")
            df1 = get_dict[s](df1, invert=invert)
            df2 = get_dict[s](df2, invert=invert)

        total = (
            (df1["x_coord"].values - df2["x_coord"].values) ** 2
            + (df1["y_coord"].values - df2["y_coord"].values) ** 2
            + (df1["z_coord"].values - df2["z_coord"].values) ** 2
        )
        return round((total.sum() / df1.shape[0]) ** 0.5, decimals)

    @staticmethod
    def _init_get_dict():
        """Initialize dictionary for filter operations."""
        return {
            "main chain": PandasPdb._get_mainchain,
            "hydrogen": PandasPdb._get_hydrogen,
            "c-alpha": PandasPdb._get_calpha,
            "carbon": PandasPdb._get_carbon,
            "heavy": PandasPdb._get_heavy,
        }

    @staticmethod
    def _read_pdb(path):
        """Read PDB file from local drive."""
        r_mode = "r"
        if path.endswith((".pdb", ".ent")):
            openf = open
        elif path.endswith(("pdb.gz", ".ent.gz")):
            r_mode = "rb"
            openf = gzip.open
        else:
            allowed_formats = ", ".join((".pdb", ".pdb.gz", ".ent", ".ent.gz"))
            raise ValueError(
                f"Wrong file format; allowed file formats are {allowed_formats}"
            )

        with openf(path, r_mode) as f:
            txt = f.read()

        if path.endswith(".gz"):
            txt = (
                txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.encode("ascii")
            )
        return path, txt

    @staticmethod
    def _fetch_pdb(pdb_code):
        """Load PDB file from rcsb.org."""
        txt = None
        url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
        try:
            response = urlopen(url)
            txt = response.read()
            txt = (
                txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.encode("ascii")
            )
        except HTTPError as e:
            print(f"HTTP Error {e.code}")
        except URLError as e:
            print(f"URL Error {e.args}")
        return url, txt

    @staticmethod
    def _fetch_af2(uniprot_id: str, af2_version: int = 3):
        """Load PDB file from https://alphafold.ebi.ac.uk/."""
        txt = None
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id.upper()}-F1-model_v{af2_version}.pdb"
        try:
            response = urlopen(url)
            txt = response.read()
            txt = (
                txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.encode("ascii")
            )
        except HTTPError as e:
            print(f"HTTP Error {e.code}")
        except URLError as e:
            print(f"URL Error {e.args}")
        return url, txt

    def _parse_header_code(self):
        """Extract header information and PDB code."""
        code, header = "", ""
        if "OTHERS" in self.df:

            header = self.df["OTHERS"][self.df["OTHERS"]["record_name"] == "HEADER"]
            if not header.empty:
                header = header["entry"].values[0]
                s = header.split()
                if s:
                    code = s[-1].lower()
        return header, code

    @staticmethod
    def _get_mainchain(df, invert):
        """Return only main chain atom entries from a DataFrame"""
        if invert:
            mc = df[
                (df["atom_name"] != "C")
                & (df["atom_name"] != "O")
                & (df["atom_name"] != "N")
                & (df["atom_name"] != "CA")
            ]
        else:
            mc = df[
                (df["atom_name"] == "C")
                | (df["atom_name"] == "O")
                | (df["atom_name"] == "N")
                | (df["atom_name"] == "CA")
            ]
        return mc

    @staticmethod
    def _get_hydrogen(df, invert):
        """Return only hydrogen atom entries from a DataFrame"""
        if invert:
            return df[(df["element_symbol"] != "H")]
        else:
            return df[(df["element_symbol"] == "H")]

    @staticmethod
    def _get_heavy(df, invert):
        """Return only heavy atom entries from a DataFrame"""
        if invert:
            return df[df["element_symbol"] == "H"]
        else:
            return df[df["element_symbol"] != "H"]

    @staticmethod
    def _get_calpha(df, invert):
        """Return c-alpha atom entries from a DataFrame"""
        if invert:
            return df[df["atom_name"] != "CA"]
        else:
            return df[df["atom_name"] == "CA"]

    @staticmethod
    def _get_carbon(df, invert):
        """Return carbon atom entries from a DataFrame"""
        if invert:
            return df[df["element_symbol"] != "C"]
        else:
            return df[df["element_symbol"] == "C"]

    @staticmethod
    def _construct_df(pdb_lines):
        """Construct DataFrames from list of PDB lines."""
        valids = tuple(pdb_records.keys())
        line_lists = {r: [] for r in valids}
        line_lists["OTHERS"] = []
        for line_num, line in enumerate(pdb_lines):
            if line.strip():
                if line.startswith(valids):
                    record = line[:6].rstrip()
                    line_ele = ["" for _ in range(len(pdb_records[record]) + 1)]
                    for idx, ele in enumerate(pdb_records[record]):
                        line_ele[idx] = line[ele["line"][0] : ele["line"][1]].strip()
                    line_ele[-1] = line_num
                    line_lists[record].append(line_ele)
                else:
                    line_lists["OTHERS"].append(
                        [line[:6].rstrip(), line[6:-1].rstrip(), line_num]
                    )

        dfs = {}
        for r in line_lists.items():
            df = pd.DataFrame(
                r[1], columns=[c["id"] for c in pdb_records[r[0]]] + ["line_idx"]
            )
            for c in pdb_records[r[0]]:
                try:
                    df[c["id"]] = df[c["id"]].astype(c["type"])
                except ValueError:
                    # expect ValueError if float/int columns are empty strings
                    df[c["id"]] = pd.Series(np.nan, index=df.index)

            dfs[r[0]] = df

        # issue a warning if no atoms have been loaded
        if (len(dfs["ATOM"]) + len(dfs["HETATM"])) == 0:
            warnings.warn(
                "No ATOM/HETATM entries have been loaded. "
                "Is the input file/text in the pdb format?"
            )

        return dfs

    def amino3to1(self, record="ATOM", residue_col="residue_name", fillna="?"):
        """Creates 1-letter amino acid codes from DataFrame

        Non-canonical amino-acids are converted as follows:
        ASH (protonated ASP) => D
        CYX (disulfide-bonded CYS) => C
        GLH (protonated GLU) => E
        HID/HIE/HIP (different protonation states of HIS) = H
        HYP (hydroxyproline) => P
        MSE (selenomethionine) => M

        Parameters
        ----------
        record : str, default: 'ATOM'
            Specifies the record DataFrame.
        residue_col : str,  default: 'residue_name'
            Column in `record` DataFrame to look for 3-letter amino acid
            codes for the conversion.
        fillna : str, default: '?'
            Placeholder string to use for unknown amino acids.

        Returns
        ---------
        pandas.DataFrame : Pandas DataFrame object consisting of two columns,
            `'chain_id'` and `'residue_name'`, where the former contains
            the chain ID of the amino acid and the latter
            contains the 1-letter amino acid code, respectively.

        """
        tmp = self.df[record]
        cmp = "placeholder"
        indices = []

        residue_number_insertion = tmp["residue_number"].astype(str) + tmp["insertion"]

        for num, ind in zip(residue_number_insertion, np.arange(tmp.shape[0])):
            if num != cmp:
                indices.append(ind)
            cmp = num

        transl = tmp.iloc[indices][residue_col].map(amino3to1dict).fillna(fillna)

        return pd.concat((tmp.iloc[indices]["chain_id"], transl), axis=1)

    def distance(self, xyz=(0.00, 0.00, 0.00), records=("ATOM", "HETATM")):
        """Computes Euclidean distance between atoms and a 3D point.

        Parameters
        ----------
        xyz : tuple, default: (0.00, 0.00, 0.00)
            X, Y, and Z coordinate of the reference center for the distance
            computation.
        records : iterable, default: ('ATOM', 'HETATM')
            Specify which record sections to consider. For example, to consider
            both protein and ligand atoms, set `records=('ATOM', 'HETATM')`.
            This setting is ignored if `df` is not set to None.
            For downward compatibility, a string argument is still supported
            but deprecated and will be removed in future versions.

        Returns
        ---------
        pandas.Series : Pandas Series object containing the Euclidean
            distance between the atoms in the record section and `xyz`.

        """

        if isinstance(records, str):
            warnings.warn(
                "Using a string as `records` argument is "
                "deprecated and will not be supported in future"
                " versions. Please use a tuple or"
                " other iterable instead",
                DeprecationWarning,
            )
            records = (records,)

        df = pd.concat(objs=[self.df[i] for i in records])

        return np.sqrt(
            np.sum(
                df[["x_coord", "y_coord", "z_coord"]].subtract(xyz, axis=1) ** 2, axis=1
            )
        )

    @staticmethod
    def distance_df(df, xyz=(0.00, 0.00, 0.00)):
        """Computes Euclidean distance between atoms and a 3D point.

        Parameters
        ----------
        df : DataFrame
            DataFrame containing entries in the `PandasPdb.df['ATOM']`
            or `PandasPdb.df['HETATM']` format for the
            the distance computation to the `xyz` reference coordinates.
        xyz : tuple, default: (0.00, 0.00, 0.00)
            X, Y, and Z coordinate of the reference center for the distance
            computation.

        Returns
        ---------
        pandas.Series : Pandas Series object containing the Euclidean
            distance between the atoms in the record section and `xyz`.

        """
        return np.sqrt(
            np.sum(
                df[["x_coord", "y_coord", "z_coord"]].subtract(xyz, axis=1) ** 2, axis=1
            )
        )

    def to_pdb(self, path, records=None, gz=False, append_newline=True):
        """Write record DataFrames to a PDB file or gzipped PDB file.

        Parameters
        ----------
        path : str
            A valid output path for the pdb file

        records : iterable, default: None
            A list of PDB record sections in
            {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'} that are to be written.
            Writes all lines to PDB if `records=None`.

        gz : bool, default: False
            Writes a gzipped PDB file if True.

        append_newline : bool, default: True
            Appends a new line at the end of the PDB file if True

        """
        if gz:
            openf = gzip.open
            w_mode = "wt"
        else:
            openf = open
            w_mode = "w"
        if not records:
            records = self.df.keys()

        dfs = {r: self.df[r].copy() for r in records if not self.df[r].empty}

        for r in dfs:
            for col in pdb_records[r]:
                dfs[r][col["id"]] = dfs[r][col["id"]].apply(col["strf"])
                dfs[r]["OUT"] = pd.Series("", index=dfs[r].index)

            for c in dfs[r].columns:
                # fix issue where coordinates with four or more digits would
                # cause issues because the columns become too wide
                if c in {"x_coord", "y_coord", "z_coord"}:
                    for idx in range(dfs[r][c].values.shape[0]):
                        if len(dfs[r][c].values[idx]) > 8:
                            dfs[r][c].values[idx] = str(dfs[r][c].values[idx]).strip()
                if c in {"line_idx", "OUT"}:
                    pass
                elif r in {"ATOM", "HETATM"} and c not in pdb_df_columns:
                    warn(
                        "Column %s is not an expected column and"
                        " will be skipped." % c
                    )
                else:
                    dfs[r]["OUT"] = dfs[r]["OUT"] + dfs[r][c]

        if pd_version < LooseVersion("0.17.0"):
            warn(
                "You are using an old pandas version (< 0.17)"
                " that relies on the old sorting syntax."
                " Please consider updating your pandas"
                " installation to a more recent version.",
                DeprecationWarning,
            )
            dfs.sort(columns="line_idx", inplace=True)

        elif pd_version < LooseVersion("0.23.0"):
            df = pd.concat(dfs)

        else:
            df = pd.concat(dfs, sort=False)

        df.sort_values(by="line_idx", inplace=True)

        with openf(path, w_mode) as f:

            s = df["OUT"].tolist()
            for idx in range(len(s)):
                if len(s[idx]) < 80:
                    s[idx] = f"{s[idx]}{' ' * (80 - len(s[idx]))}"
            to_write = "\n".join(s)
            f.write(to_write)
            if append_newline:
                f.write("\n")

    def parse_sse(self):
        """Parse secondary structure elements"""
        raise NotImplementedError

    def get_model_start_end(self) -> pd.DataFrame:
        """Get the start and end of the models contained in the PDB file.

        Extracts model start and end line indexes based
          on lines labelled 'OTHERS' during parsing.

        Returns
        ---------
        pandas.DataFrame : Pandas DataFrame object containing
          the start and end line indexes of the models.
        """

        other_records = self.df["OTHERS"]

        idxs = other_records.loc[other_records["record_name"] == "MODEL"].copy()
        ends = other_records.loc[other_records["record_name"] == "ENDMDL"]
        idxs.columns = ["record_name", "model_idx", "start_idx"]
        idxs.loc[:, "end_idx"] = ends.line_idx.values
        # If structure only contains 1 model, create a dummy df mapping all lines to model_idx 1
        if len(idxs) == 0:
            n_lines = len(self.pdb_text.splitlines())
            idxs = pd.DataFrame(
                [
                    {
                        "record_name": "MODEL",
                        "model_idx": 1,
                        "start_idx": 0,
                        "end_idx": n_lines,
                    }
                ]
            )

        return idxs

    def label_models(self):
        """Adds a column (`"model_id"`) to the underlying
        DataFrames containing the model number."""
        idxs = self.get_model_start_end()
        # Label ATOMS
        if "ATOM" in self.df.keys():
            pdb_df = self.df["ATOM"]
            idx_map = np.piecewise(
                np.zeros(len(pdb_df)),
                [
                    (pdb_df.line_idx.values >= start_idx)
                    & (pdb_df.line_idx.values <= end_idx)
                    for start_idx, end_idx in zip(
                        idxs.start_idx.values, idxs.end_idx.values
                    )
                ],
                idxs.model_idx,
            )
            self.df["ATOM"]["model_id"] = idx_map.astype(int)
        # LABEL HETATMS
        if "HETATM" in self.df.keys():
            pdb_df = self.df["HETATM"]
            idx_map = np.piecewise(
                np.zeros(len(pdb_df)),
                [
                    (pdb_df.line_idx.values >= start_idx)
                    & (pdb_df.line_idx.values <= end_idx)
                    for start_idx, end_idx in zip(
                        idxs.start_idx.values, idxs.end_idx.values
                    )
                ],
                idxs.model_idx,
            )
            self.df["HETATM"]["model_id"] = idx_map.astype(int)
        if "ANISOU" in self.df.keys():
            pdb_df = self.df["ANISOU"]
            idx_map = np.piecewise(
                np.zeros(len(pdb_df)),
                [
                    (pdb_df.line_idx.values >= start_idx)
                    & (pdb_df.line_idx.values <= end_idx)
                    for start_idx, end_idx in zip(
                        idxs.start_idx.values, idxs.end_idx.values
                    )
                ],
                idxs.model_idx,
            )
            self.df["ANISOU"]["model_id"] = idx_map.astype(int)
        return self

    def get_model(self, model_index: int) -> PandasPdb:
        """Returns a new PandasPDB object with the dataframes subset to the given model index.

        Parameters
        ----------
        model_index : int
            An integer representing the model index to subset to.

        Returns
        ---------
        pandas_pdb.PandasPdb : A new PandasPdb object containing the
          structure subsetted to the given model.
        """

        biopandas_structure = deepcopy(self)
        biopandas_structure.label_models()

        if "ATOM" in biopandas_structure.df.keys():
            biopandas_structure.df["ATOM"] = biopandas_structure.df["ATOM"].loc[
                biopandas_structure.df["ATOM"]["model_id"] == model_index
            ]
        if "HETATM" in biopandas_structure.df.keys():
            biopandas_structure.df["HETATM"] = biopandas_structure.df["HETATM"].loc[
                biopandas_structure.df["HETATM"]["model_id"] == model_index
            ]
        if "ANISOU" in biopandas_structure.df.keys():
            biopandas_structure.df["ANISOU"] = biopandas_structure.df["ANISOU"].loc[
                biopandas_structure.df["ANISOU"]["model_id"] == model_index
            ]
        return biopandas_structure

    def get_models(self, model_indices: List[int]) -> PandasPdb:
        """Returns a new PandasPDB object with the dataframes subset to the given model index.

        Parameters
        ----------
        model_indices : List[int]
            A list representing the model indexes to subset to.

        Returns
        ---------
        pandas_pdb.PandasPdb : A new PandasPdb object
          containing the structure subsetted to the given model.
        """

        biopandas_structure = deepcopy(self)
        biopandas_structure.label_models()

        if "ATOM" in biopandas_structure.df.keys():
            biopandas_structure.df["ATOM"] = biopandas_structure.df["ATOM"].loc[
                [
                    x in model_indices
                    for x in biopandas_structure.df["ATOM"]["model_id"].tolist()
                ]
            ]
        if "HETATM" in biopandas_structure.df.keys():
            biopandas_structure.df["HETATM"] = biopandas_structure.df["HETATM"].loc[
                [
                    x in model_indices
                    for x in biopandas_structure.df["HETATM"]["model_id"].tolist()
                ]
            ]
        if "ANISOU" in biopandas_structure.df.keys():
            biopandas_structure.df["ANISOU"] = biopandas_structure.df["ANISOU"].loc[
                [
                    x in model_indices
                    for x in biopandas_structure.df["ANISOU"]["model_id"].tolist()
                ]
            ]
        return biopandas_structure

    def to_pdb_stream(self, records: tuple[str] = ("ATOM", "HETATM")) -> StringIO:
        """Writes a PDB dataframe to a stream.

        Parameters
        ------------
        records : iterable, default: ('ATOM', 'HETATM')
            Iterable of record names to save to stream. Any of `["ATOM", "HETATM", "OTHERS"]`.

        Returns
        --------
        io.StringIO : Filestream of PDB file.
        """

        df = self.df.copy()
        df = pd.concat([df[a] for a in records])
        if "model_id" in df.columns:
            df = df.drop(columns=["model_id"])
        df["residue_number"] = pd.to_numeric(df.residue_number, errors="coerce")
        records = [r.strip() for r in list(set(df.record_name))]
        dfs = {r: df.loc[df.record_name == r] for r in records}

        for r in dfs:
            for col in pdb_records[r]:
                dfs[r][col["id"]] = dfs[r][col["id"]].apply(col["strf"])
                dfs[r]["OUT"] = pd.Series("", index=dfs[r].index)

            for c in dfs[r].columns:
                # fix issue where coordinates with four or more digits would
                # cause issues because the columns become too wide
                if c in {"x_coord", "y_coord", "z_coord"}:
                    for idx in range(dfs[r][c].values.shape[0]):
                        if len(dfs[r][c].values[idx]) > 8:
                            dfs[r][c].values[idx] = str(dfs[r][c].values[idx]).strip()

                if c not in {"line_idx", "OUT"}:
                    dfs[r]["OUT"] = dfs[r]["OUT"] + dfs[r][c]

        df = pd.concat(dfs, sort=False)
        df.sort_values(by="line_idx", inplace=True)

        output = StringIO()
        s = df["OUT"].tolist()
        for idx in range(len(s)):
            if len(s[idx]) < 80:
                s[idx] = f"{s[idx]}{' ' * (80 - len(s[idx]))}"
        to_write = "\n".join(s)
        output.write(to_write)
        output.write("\n")
        output.seek(0)
        return output

    def gyradius(self, records: tuple[str] = ("ATOM",), decimals: int = 4) -> float:
        """Compute the Radius of Gyration of a molecule

        Parameters
        ----------
        records : iterable, default: ("ATOM",)
            Records from PandasPdb object for which to calculate the radius of gyration.
            Any of `("ATOM", "HETATM")`.

        decimals : int, default: 4
            Specifies the number of decimal places to round the final value to.

        Returns
        ---------
        rg : float
            Radius of Gyration of df in Angstrom

        """
        if isinstance(records, str):
            warnings.warn(
                "Using a string as `records` argument is "
                "deprecated and will not be supported in future "
                "versions. Please use a tuple or "
                "other iterable instead",
                DeprecationWarning,
            )
            records = (records,)

        if len(records) > 1:
            df = pd.concat(
                objs=[
                    self.df[record][["x_coord", "y_coord", "z_coord", "element_symbol"]]
                    for record in records
                ]
            )
        else:
            df = self.df[records[0]]

        coords = df[["x_coord", "y_coord", "z_coord"]].to_numpy()
        masses = (
            df["element_symbol"].map(lambda atom: ATOMIC_MASSES.get(atom, 0)).to_numpy()
        )
        total_mass = masses.sum()
        center_of_mass = (masses[:, None] * coords).sum(axis=0) / total_mass
        distances = np.linalg.norm(coords - center_of_mass, axis=1)
        rg = np.sqrt((distances**2 * masses).sum() / total_mass)
        return round(rg, decimals)
