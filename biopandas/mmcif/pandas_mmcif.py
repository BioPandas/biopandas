"""Class for working with MMCIF files."""

# BioPandas
# Authors: Arian Jamasb <arian@jamasb.io>,
# Authors: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas
from __future__ import annotations
import gzip
import sys
import copy
import warnings
from typing import Dict, List, Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd
from looseversion import LooseVersion

from ..pdb.engines import amino3to1dict
from ..pdb.pandas_pdb import PandasPdb
from .engines import (ANISOU_DF_COLUMNS, MMCIF_PDB_COLUMN_MAP,
                      MMCIF_PDB_NONEFIELDS, PDB_COLUMN_ORDER, mmcif_col_types)
from .mmcif_parser import load_cif_data

pd_version = LooseVersion(pd.__version__)


class PandasMmcif:
    def __init__(self, use_auth: bool = True):
        self._df = None
        self.mmcif_text = ""
        self.header = ""
        self.code = ""
        self.mmcif_path = ""
        self.auth = use_auth
        self._get_dict = {}

    @property
    def df(self):
        """Acccess dictionary of pandas DataFrames for PDB record sections."""
        return self._df

    @df.setter
    def df(self, value):
        """Assign a new value to the pandas DataFrame"""
        raise AttributeError(
            "Please use `PandasMmcif._df = ... ` instead\n"
            "of `PandasMmcif.df = ... ` if you are sure that\n"
            "you want to overwrite the `df` attribute."
        )

    def read_mmcif(self, path):
        """Read MMCIF files (unzipped or gzipped) from local drive

        Attributes
        ----------
        path : Union[str, os.PathLike]
            Path to the MMCIF file in .cif format or gzipped format (.cif.gz).

        Returns
        ---------
        self

        """
        self.mmcif_path, self.pdb_text = self._read_mmcif(path=str(path))
        self._df = self._construct_df(text=self.pdb_text)
        # self.header, self.code = self._parse_header_code() #TODO: implement
        self.code = self.data["entry"]["id"][0].lower()
        return self
    
    def label_models(self):
        """Adds a column ("model_id") to the underlying
        DataFrames containing the model number."""
        if "ATOM" in self.df.keys():
            self.df["ATOM"]["model_id"] = self.df["ATOM"]["pdbx_PDB_model_num"]
        if "HETATM" in self.df.keys():
            self.df["HETATM"]["model_id"] = self.df["HETATM"]["pdbx_PDB_model_num"]
        return self
    
    def get_model(self, model_index: int) -> PandasMmcif:
        """Returns a new PandasMmcif object with the dataframes subset to the
        given model index.

        Parameters
        ----------
        model_index : int
            An integer representing the model index to subset to.

        Returns
        ---------
        pandas_pdb.PandasPdb : A new PandasMMcif object containing the
            structure subsetted to the given model.
        """

        biopandas_structure = copy.deepcopy(self)
        if "ATOM" in biopandas_structure.df.keys():
            biopandas_structure.df["ATOM"] = biopandas_structure.df["ATOM"].loc[biopandas_structure.df["ATOM"]["pdbx_PDB_model_num"] == model_index]
        if "HETATM" in biopandas_structure.df.keys():
            biopandas_structure.df["HETATM"] = biopandas_structure.df["HETATM"].loc[
                biopandas_structure.df["HETATM"]["pdbx_PDB_model_num"] == model_index
            ]
        return biopandas_structure

    def get_models(self, model_indices: List[int]) -> PandasMmcif:
        """Returns a new PandasMmcif object with the dataframes subset to the
        given model index.

        Parameters
        ----------
        model_indices : List[int]
            A list representing the model indexes to subset to.

        Returns
        ---------
        pandas_pdb.PandasMmtf : A new PandasMmcif object
            containing the structure subsetted to the given model.
        """

        biopandas_structure = copy.deepcopy(self)

        if "ATOM" in biopandas_structure.df.keys():
            biopandas_structure.df["ATOM"] = biopandas_structure.df["ATOM"].loc[
                [x in model_indices for x in biopandas_structure.df["ATOM"]["pdbx_PDB_model_num"].tolist()]
            ]
        if "HETATM" in biopandas_structure.df.keys():
            biopandas_structure.df["HETATM"] = biopandas_structure.df["HETATM"].loc[
                [x in model_indices for x in biopandas_structure.df["HETATM"]["pdbx_PDB_model_num"].tolist()]
            ]
        return biopandas_structure

    def fetch_mmcif(
        self,
        pdb_code: Optional[str] = None,
        uniprot_id: Optional[str] = None,
        source: str = "pdb",
    ):
        """Fetches mmCIF file contents from the Protein Databank at rcsb.org or AlphaFold database at https://alphafold.ebi.ac.uk/.
        .

                Parameters
                ----------
                pdb_code : str, optional
                    A 4-letter PDB code, e.g., `"3eiy"` to retrieve structures from the PDB. Defaults to `None`.

                uniprot_id : str, optional
                    A UniProt Identifier, e.g., `"Q5VSL9"` to retrieve structures from the AF2 database. Defaults to `None`.

                source : str
                    The source to retrieve the structure from
                    (`"pdb"`, `"alphafold2-v3"` or `"alphafold2-v4"`). Defaults to `"pdb"`.

                Returns
                ---------
                self

        """
        # Sanitize input
        invalid_input_identifier_1 = pdb_code is None and uniprot_id is None
        invalid_input_identifier_2 = (
            pdb_code is not None and uniprot_id is not None
        )
        invalid_input_combination_1 = (
            uniprot_id is not None and source == "pdb"
        )
        invalid_input_combination_2 = pdb_code is not None and source in {
            "alphafold2-v3",
            "alphafold2-v4",
        }

        if invalid_input_identifier_1 or invalid_input_identifier_2:
            raise ValueError(
                "Please provide either a PDB code or a UniProt ID."
            )

        if invalid_input_combination_1:
            raise ValueError(
                "Please use a 'pdb_code' instead of 'uniprot_id' for source='pdb'."
            )
        elif invalid_input_combination_2:
            raise ValueError(
                f"Please use a 'uniprot_id' instead of 'pdb_code' for source={source}."
            )

        if source == "pdb":
            self.mmcif_path, self.mmcif_text = self._fetch_mmcif(pdb_code)
        elif source == "alphafold2-v3":
            af2_version = 3
            self.mmcif_path, self.mmcif_text = self._fetch_af2(
                uniprot_id, af2_version
            )
        elif source == "alphafold2-v4":
            af2_version = 4
            self.mmcif_path, self.mmcif_text = self._fetch_af2(
                uniprot_id, af2_version
            )
        else:
            raise ValueError(
                f"Invalid source: {source}."
                " Please use one of 'pdb', 'alphafold2-v3' or 'alphafold2-v4'."
            )

        self._df = self._construct_df(text=self.mmcif_text)
        return self

    def _construct_df(self, text: str):
        data = load_cif_data(text)
        data = data[list(data.keys())[0]]
        self.data = data
        df: Dict[str, pd.DataFrame] = {}
        full_df = pd.DataFrame.from_dict(
            data["atom_site"], orient="index"
        ).transpose()
        full_df = full_df.astype(mmcif_col_types, errors="ignore")
        df["ATOM"] = pd.DataFrame(full_df[full_df.group_PDB == "ATOM"])
        df["HETATM"] = pd.DataFrame(full_df[full_df.group_PDB == "HETATM"])
        try:
            df["ANISOU"] = pd.DataFrame(data["atom_site_anisotrop"])
        except KeyError:
            df["ANISOU"] = pd.DataFrame(columns=ANISOU_DF_COLUMNS)
        return df

    @staticmethod
    def _fetch_mmcif(pdb_code):
        """Load MMCIF file from rcsb.org."""
        txt = None
        url = f"https://files.rcsb.org/download/{pdb_code.lower()}.cif"
        try:
            response = urlopen(url)
            txt = response.read()
            txt = (
                txt.decode("utf-8")
                if sys.version_info[0] >= 3
                else txt.encode("ascii")
            )
        except HTTPError as e:
            print(f"HTTP Error {e.code}")
        except URLError as e:
            print(f"URL Error {e.args}")
        return url, txt

    @staticmethod
    def _fetch_af2(uniprot_id: str, af2_version: int = 3):
        """Load MMCIF file from https://alphafold.ebi.ac.uk/."""
        txt = None
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id.upper()}-F1-model_v{af2_version}.cif"

        try:
            response = urlopen(url)
            txt = response.read()
            txt = (
                txt.decode("utf-8")
                if sys.version_info[0] >= 3
                else txt.encode("ascii")
            )
        except HTTPError as e:
            print(f"HTTP Error {e.code}")
        except URLError as e:
            print(f"URL Error {e.args}")
        return url, txt

    @staticmethod
    def _read_mmcif(path):
        """Read MMCIF file from local drive."""
        r_mode = "r"
        if path.endswith((".cif", ".mmcif")):
            openf = open
        elif path.endswith((".cif.gz", ".mmcif.gz")):
            r_mode = "rb"
            openf = gzip.open
        else:
            allowed_formats = ", ".join(
                (".cif", ".cif.gz", ".mmcif", ".mmcif.gz")
            )
            raise ValueError(
                f"Wrong file format; allowed file formats are {allowed_formats}"
            )

        with openf(path, r_mode) as f:
            txt = f.read()

        if path.endswith(".gz"):
            txt = (
                txt.decode("utf-8")
                if sys.version_info[0] >= 3
                else txt.encode("ascii")
            )
        return path, txt

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

    @staticmethod
    def _get_mainchain(
        df: pd.DataFrame, invert: bool = False, atom_col: str = "auth_atom_id"
    ) -> pd.DataFrame:
        """Return only main chain atom entries from a DataFrame"""
        return (
            df[
                (df[atom_col] != "C")
                & (df[atom_col] != "O")
                & (df[atom_col] != "N")
                & (df[atom_col] != "CA")
            ]
            if invert
            else df[
                (df[atom_col] == "C")
                | (df[atom_col] == "O")
                | (df[atom_col] == "N")
                | (df[atom_col] == "CA")
            ]
        )

    @staticmethod
    def _get_hydrogen(df, invert):
        """Return only hydrogen atom entries from a DataFrame"""
        return (
            df[(df["type_symbol"] != "H")]
            if invert
            else df[(df["type_symbol"] == "H")]
        )

    @staticmethod
    def _get_heavy(df, invert):
        """Return only heavy atom entries from a DataFrame"""
        return (
            df[df["type_symbol"] == "H"]
            if invert
            else df[df["type_symbol"] != "H"]
        )

    @staticmethod
    def _get_calpha(df, invert, atom_col: str = "auth_atom_id"):
        """Return c-alpha atom entries from a DataFrame"""
        return df[df[atom_col] != "CA"] if invert else df[df[atom_col] == "CA"]

    @staticmethod
    def _get_carbon(df, invert):
        """Return carbon atom entries from a DataFrame"""
        return (
            df[df["type_symbol"] != "C"]
            if invert
            else df[df["type_symbol"] == "C"]
        )

    def amino3to1(
        self,
        record: str = "ATOM",
        residue_col: str = "auth_comp_id",
        residue_number_col: str = "auth_seq_id",
        chain_col: str = "auth_asym_id",
        fillna: str = "?",
    ):
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

        residue_number_insertion = (
            tmp[residue_number_col].astype(str) + tmp["pdbx_PDB_ins_code"]
        )

        for num, ind in zip(residue_number_insertion, np.arange(tmp.shape[0])):
            if num != cmp:
                indices.append(ind)
            cmp = num

        transl = (
            tmp.iloc[indices][residue_col].map(amino3to1dict).fillna(fillna)
        )

        return pd.concat((tmp.iloc[indices][chain_col], transl), axis=1)

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
        get_dict = PandasMmcif._init_get_dict()
        if s:
            if s not in get_dict.keys():
                raise AttributeError(f"s must be in {get_dict.keys()} or None")
            df1 = get_dict[s](df1, invert=invert)
            df2 = get_dict[s](df2, invert=invert)

        total = (
            (df1["Cartn_x"].values - df2["Cartn_x"].values) ** 2
            + (df1["Cartn_y"].values - df2["Cartn_y"].values) ** 2
            + (df1["Cartn_z"].values - df2["Cartn_z"].values) ** 2
        )
        return round((total.sum() / df1.shape[0]) ** 0.5, 4)

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
                df[["Cartn_x", "Cartn_y", "Cartn_z"]].subtract(xyz, axis=1)
                ** 2,
                axis=1,
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
                df[["Cartn_x", "Cartn_y", "Cartn_z"]].subtract(xyz, axis=1)
                ** 2,
                axis=1,
            )
        )

    @staticmethod
    def _init_get_dict():
        """Initialize dictionary for filter operations."""
        return {
            "main chain": PandasMmcif._get_mainchain,
            "hydrogen": PandasMmcif._get_hydrogen,
            "c-alpha": PandasMmcif._get_calpha,
            "carbon": PandasMmcif._get_carbon,
            "heavy": PandasMmcif._get_heavy,
        }

    def read_mmcif_from_list(self, mmcif_lines):
        """Reads mmCIF file from a list into DataFrames

        Attributes
        ----------
        pdb_lines : list
            A list of lines containing the mmCIF file contents.

        Returns
        ---------
        self

        """
        self.pdb_text = "".join(mmcif_lines)
        self._df = self._construct_df(mmcif_lines)
        # self.header, self.code = self._parse_header_code()
        self.code = self.data["entry"]["id"][0].lower()
        return self

    def convert_to_pandas_pdb(
        self,
        offset_chains: bool = True,
        records: List[str] = ["ATOM", "HETATM"],
    ) -> PandasPdb:
        """Returns a PandasPdb object with the same data as the PandasMmcif
        object.

        Attributes
        ----------
        offset_chains: bool
            Whether or not to offset atom numbering based on number of chains.
            This can arise due to the presence of TER records in PDBs which are
            not found in mmCIFs.
        records: List[str]
            List of record types to save. Any of ["ATOM", "HETATM", "OTHERS"].
            Defaults to ["ATOM", "HETATM"].

        """
        pandaspdb = PandasPdb()

        for a in records:
            try:
                dfa = self.df[a]
                # keep only those fields found in pdb
                dfa = dfa[MMCIF_PDB_COLUMN_MAP.keys()]
                # rename fields
                dfa = dfa.rename(columns=MMCIF_PDB_COLUMN_MAP)
                # add empty fields
                for i in MMCIF_PDB_NONEFIELDS:
                    dfa[i] = ""
                dfa["charge"] = np.nan
                # reorder columns to PandasPdb order
                dfa = dfa[PDB_COLUMN_ORDER]
                pandaspdb.df[a] = dfa
            except KeyError:  # Some entries may not have an ANISOU
                continue

        # update line_idx
        pandaspdb.df["ATOM"]["line_idx"] = pandaspdb.df["ATOM"].index.values
        pandaspdb.df["HETATM"]["line_idx"] = pandaspdb.df["HETATM"].index

        # Update atom numbers
        if offset_chains:
            offsets = (
                pandaspdb.df["ATOM"]["chain_id"].astype("category").cat.codes
            )
            pandaspdb.df["ATOM"]["atom_number"] = (
                pandaspdb.df["ATOM"]["atom_number"] + offsets
            )
            hetatom_offset = offsets.max() + 1
            pandaspdb.df["HETATM"]["atom_number"] = (
                pandaspdb.df["HETATM"]["atom_number"] + hetatom_offset
            )

        return pandaspdb