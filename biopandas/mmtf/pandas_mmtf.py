"""Class for working with MMTF files."""

import gzip
import warnings
from distutils.version import LooseVersion
from typing import Any, Dict, List, Union
from warnings import warn

import numpy as np
import pandas as pd

from mmtf import MMTFDecoder, fetch, parse

from ..pdb.engines import amino3to1dict, pdb_df_columns, pdb_records

pd_version = LooseVersion(pd.__version__)


class PandasMmtf:
    def __init__(self):
        self._df: Dict[str, pd.DataFrame] = {}
        self.mmtf = ""
        self.header = ""
        self.code = ""
        self._get_dict = {}
        self.mmtf_path = ""

    @property
    def df(self):
        """Access dictionary of pandas DataFrames for MMTF record sections."""
        return self._df

    @df.setter
    def df(self, value: Any):
        """Assign a new value to the pandas DataFrame"""
        raise AttributeError(
            "Please use `PandasMmtf._df = ... ` instead\n"
            "of `PandasMmtf.df = ... ` if you are sure that\n"
            "you want to overwrite the `df` attribute."
        )
        # self._df = value

    def read_mmtf(self, filename: str):
        self.mmtf = parse(filename)
        self.mmtf_path = filename
        df = self._mmtf_to_df(self.mmtf)
        self._df["ATOM"] = df.loc[df.record_name == "ATOM"]
        self._df["HETATM"] = df.loc[df.record_name == "HETATM"]

    def fetch_mmtf(self, pdb_code: str):
        self.code = pdb_code
        self.mmtf = fetch(pdb_code)
        df = self._mmtf_to_df(self.mmtf)
        self._df["ATOM"] = df.loc[df.record_name == "ATOM"]
        self._df["HETATM"] = df.loc[df.record_name == "HETATM"]

    @staticmethod
    def _mmtf_to_df(mmtf_obj: MMTFDecoder) -> pd.DataFrame:
        return mmtf_to_df(mmtf_obj)

    def impute_element(self, records=("ATOM", "HETATM"), inplace=False):
        """Impute element_symbol from atom_name section.

        Parameters
        ----------
        records : iterable, default: ('ATOM', 'HETATM')
            Coordinate sections for which the element symbols should be
            imputed.

        inplace : bool, (default: False
            Performs the operation in-place if True and returns a copy of the
            MMTF DataFrame otherwise.

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

    @staticmethod
    def rmsd(df1, df2, s=None, invert=False):
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

        Returns
        ---------
        rmsd : float
            Root Mean Square Deviation between df1 and df2

        """
        if df1.shape[0] != df2.shape[0]:
            raise AttributeError("DataFrames have unequal lengths")
        get_dict = PandasMmtf._init_get_dict()
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
        return round((total.sum() / df1.shape[0]) ** 0.5, 4)

    @staticmethod
    def _init_get_dict():
        """Initialize dictionary for filter operations."""
        return {
            "main chain": PandasMmtf._get_mainchain,
            "hydrogen": PandasMmtf._get_hydrogen,
            "c-alpha": PandasMmtf._get_calpha,
            "carbon": PandasMmtf._get_carbon,
            "heavy": PandasMmtf._get_heavy,
        }

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
            Specfies the record DataFrame.
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
            DataFrame containing entries in the `PandasMmtf.df['ATOM']`
            or `PandasMmtf.df['HETATM']` format for the
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


def fetch_mmtf(pdb_code: str) -> pd.DataFrame:
    """Returns a dataframe from a PDB code.

    :param pdb_code: 4-letter PDB accession code
    :type pdb_code: str
    :return: Dataframe of protein structure
    :rtype: pd.DataFrame
    """
    df = fetch(pdb_code)
    return mmtf_to_df(df)


def parse_mmtf(file_path: str) -> pd.DataFrame:
    """Parses an MMTF file from disk and returns a pandas DataFrame.

    :param file_path: Path to mmtf file.
    :type file_path: str
    :return: Dataframe of protein structure.
    :rtype: pd.DataFrame
    """
    df = parse(file_path)
    return mmtf_to_df(df)


def mmtf_to_df(mmtf_obj: MMTFDecoder) -> pd.DataFrame:
    data: Dict[str, Union[List[str], List[int], List[float], np.ndarray]] = {
        "record_name": [],
        "residue_name": [],
        "atom_name": [],
        "element_symbol": [],
        "charge": [],
        "x_coord": mmtf_obj.x_coord_list,
        "y_coord": mmtf_obj.y_coord_list,
        "z_coord": mmtf_obj.z_coord_list,
        "alt_loc": mmtf_obj.alt_loc_list,
        "b_factor": mmtf_obj.b_factor_list,
        "insertion": [],
        "residue_number": [],
        "occupancy": mmtf_obj.occupancy_list,
        "chain_id": [],
    }

    chain_indices = mmtf_obj.groups_per_chain

    for i, _ in enumerate(chain_indices):
        if i == 0:
            continue
        else:
            chain_indices[i] = chain_indices[i] + chain_indices[i - 1]

    ch_idx = 0
    for idx, i in enumerate(mmtf_obj.group_type_list):
        res = mmtf_obj.group_list[i]
        if idx == chain_indices[ch_idx]:
            ch_idx += 1
        record = "HETATM" if res["chemCompType"] == "NON-POLYMER" else "ATOM"
        for _ in res["atomNameList"]:
            data["residue_name"].append(res["groupName"])
            data["residue_number"].append(mmtf_obj.group_id_list[i])
            data["chain_id"].append(mmtf_obj.chain_name_list[ch_idx])
            data["record_name"].append(record)
            data["insertion"].append(mmtf_obj.ins_code_list[idx])
        data["atom_name"].append(res["atomNameList"])
        data["element_symbol"].append(res["elementList"])
        data["charge"].append(res["formalChargeList"])

    for k, v in data.items():
        if k in [
            "residue_name",
            "x_coord",
            "y_coord",
            "z_coord",
            "alt_loc",
            "b_factor",
            "residue_number",
            "occupancy",
            "record_name",
            "insertion",
        ]:
            continue
        data[k] = [i for sublist in v for i in sublist]

    df = pd.DataFrame.from_dict(data)
    df["atom_number"] = df.index + 1
    return df
