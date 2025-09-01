"""Class for working with MMTF files."""

from __future__ import annotations

import os
import copy
import gzip
import warnings
from typing import Any, Dict, List, Union
from warnings import warn

import numpy as np
import pandas as pd
from looseversion import LooseVersion
from mmtf import MMTFDecoder, MMTFEncoder, fetch, parse, parse_gzip

from biopandas.constants import protein_letters_3to1_extended

from ..pdb.engines import amino3to1dict, pdb_df_columns, pdb_records

pd_version = LooseVersion(pd.__version__)


class PandasMmtf(object):
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

    def read_mmtf(self, filename: Union[str, os.PathLike]):
        filename = str(filename)
        if filename.endswith(".gz"):
            self.mmtf = parse_gzip(filename)
        else:
            self.mmtf = parse(filename)
        self.mmtf_path = filename
        df = self._mmtf_to_df(self.mmtf)
        self._df["ATOM"] = df.loc[df.record_name == "ATOM"]
        self._df["HETATM"] = df.loc[df.record_name == "HETATM"]
        return self

    def fetch_mmtf(self, pdb_code: str):
        raise DeprecationWarning("PDB No longer serves MMTF files.")

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
            t[sec]["element_symbol"] = t[sec][
                ["atom_name", "element_symbol"]
            ].apply(lambda x: x[0][1] if len(x[1]) == 3 else x[0][0], axis=1)
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

        residue_number_insertion = (
            tmp["residue_number"].astype(str) + tmp["insertion"]
        )

        for num, ind in zip(residue_number_insertion, np.arange(tmp.shape[0])):
            if num != cmp:
                indices.append(ind)
            cmp = num

        transl = (
            tmp.iloc[indices][residue_col].map(amino3to1dict).fillna(fillna)
        )

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
                df[["x_coord", "y_coord", "z_coord"]].subtract(xyz, axis=1)
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
                df[["x_coord", "y_coord", "z_coord"]].subtract(xyz, axis=1)
                ** 2,
                axis=1,
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
                            dfs[r][c].values[idx] = str(
                                dfs[r][c].values[idx]
                            ).strip()
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

    def to_mmtf(self, path, records=("ATOM", "HETATM")):
        """Write record DataFrames to an MMTF file.

        Parameters
        ----------
        path : str
            A valid output path for the mmtf file

        records: tuple(str):
            A tuple of records to write. Defaults to ("ATOM". "HETATM")

        """
        df = pd.concat(objs=[self.df[i] for i in records])
        return write_mmtf(df, path)

    def get_model(self, model_index: int) -> PandasMmtf:
        """Returns a new PandasPDB object with the dataframes subset to the
        given model index.

        Parameters
        ----------
        model_index : int
            An integer representing the model index to subset to.

        Returns
        ---------
        pandas_pdb.PandasPdb : A new PandasPdb object containing the
            structure subsetted to the given model.
        """

        biopandas_structure = copy.deepcopy(self)

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

    def get_models(self, model_indices: List[int]) -> PandasMmtf:
        """Returns a new PandasMmtf object with the dataframes subset to the
        given model index.

        Parameters
        ----------
        model_indices : List[int]
            A list representing the model indexes to subset to.

        Returns
        ---------
        pandas_pdb.PandasMmtf : A new PandasPdb object
            containing the structure subsetted to the given model.
        """

        biopandas_structure = copy.deepcopy(self)

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
    df = (
        parse_gzip(file_path)
        if file_path.endswith(".gz")
        else parse(file_path)
    )
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
        "atom_number": mmtf_obj.atom_id_list,
        "model_id": [],
    }

    chain_indices = mmtf_obj.groups_per_chain

    for i, _ in enumerate(chain_indices):
        if i == 0:
            continue
        else:
            chain_indices[i] = chain_indices[i] + chain_indices[i - 1]
    model_indices = mmtf_obj.chains_per_model
    model_indices = [
        sum(model_indices[: i + 1]) for i in range(len(model_indices))
    ]
    ch_idx = 0

    entity_types = {}
    for i in mmtf_obj.entity_list:
        for chain in i["chainIndexList"]:
            entity_types[chain] = i["type"]

    ch_idx_iter = iter(sorted(list(entity_types.keys())))
    ch_idx = next(ch_idx_iter)
    for idx, i in enumerate(mmtf_obj.group_type_list):
        res = mmtf_obj.group_list[i]
        # record = "HETATM" if res["chemCompType"] == "NON-POLYMER" else "ATOM"
        # record = (
        #    "ATOM"
        #    if res["chemCompType"] in ["L-PEPTIDE LINKING", "PEPTIDE LINKING"]
        #    else "HETATM"
        # )
        if idx == chain_indices[ch_idx]:
            # ch_idx += 1
            ch_idx = next(ch_idx_iter)
        record = "ATOM" if entity_types[ch_idx] == "polymer" else "HETATM"

        for _ in res["atomNameList"]:
            data["residue_name"].append(res["groupName"])
            data["residue_number"].append(mmtf_obj.group_id_list[idx])
            # data["chain_id"].append([mmtf_obj.chain_name_list[ch_idx]])
            data["chain_id"].append([mmtf_obj.chain_name_list[ch_idx]])
            data["model_id"].append(
                int(np.argwhere(np.array(model_indices) > ch_idx)[0]) + 1
            )
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
            "atom_number",
            "model_id",
        ]:
            continue
        data[k] = [i for sublist in v for i in sublist]

    df = pd.DataFrame.from_dict(data).sort_values(
        by=["model_id", "atom_number"]
    )
    df.alt_loc = df.alt_loc.str.replace("\x00", "")
    df.insertion = df.insertion.str.replace("\x00", "")
    return df


def _seq1(seq, charmap: Dict[str, str], undef_code="X"):
    # sourcery skip: dict-assign-update-to-union
    """Convert protein sequence from three-letter to one-letter code.
    The single required input argument 'seq' should be a protein sequence
    using three-letter codes, either as a Python string or as a Seq or
    MutableSeq object.
    This function returns the amino acid sequence as a string using the one
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters "B" for "Asx", "J" for "Xle", "X" for "Xaa", "U" for
    "Sel", and "O" for "Pyl") plus "*" for a terminator given the "Ter" code.
    Any unknown character (including possible gap characters), is changed
    into '-' by default.
    e.g.
    >>> from Bio.SeqUtils import seq1
    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer")
    'MAIVMGRWKGAR*'
    The input is case insensitive, e.g.
    >>> from Bio.SeqUtils import seq1
    >>> seq1("METalaIlEValMetGLYArgtRplysGlyAlaARGTer")
    'MAIVMGRWKGAR*'
    You can set a custom translation of the codon termination code using the
    dictionary "custom_map" argument (defaulting to {'Ter': '*'}), e.g.
    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAla***", custom_map={"***": "*"})
    'MAIVMGRWKGA*'
    You can also set a custom translation for non-amino acid characters, such
    as '-', using the "undef_code" argument, e.g.
    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer", undef_code='?')
    'MAIVMGRWKGA??R*'
    If not given, "undef_code" defaults to "X", e.g.
    >>> seq1("MetAlaIleValMetGlyArgTrpLysGlyAla------ArgTer")
    'MAIVMGRWKGAXXR*'
    """
    if charmap is None:
        charmap = {"Ter": "*"}
    # reverse map of threecode
    # upper() on all keys to enable caps-insensitive input seq handling
    onecode = {k.upper(): v for k, v in charmap.items()}
    # add the given termination codon code and custom maps
    onecode.update((k.upper(), v) for k, v in charmap.items())
    seqlist = [seq[3 * i:3 * (i + 1)] for i in range(len(seq) // 3)]
    return "".join(onecode.get(aa.upper(), undef_code) for aa in seqlist)


def write_mmtf(df: pd.DataFrame, file_path: str):
    """Writes a biopandas dataframe to an MMTF file.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to write
    file_path : str, default: '?'
        Path to output file.
    """
    count_models, count_chains, count_groups, count_atoms = 0, 0, 0, 0
    # Check if the input is a valid BioPandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise TypeError("The input must be a BioPandas DataFrame.")

    # Initialize MMTF encoder
    encoder = MMTFEncoder()
    encoder.init_structure(
        total_num_bonds=0,
        total_num_atoms=0,
        total_num_groups=0,
        total_num_chains=0,
        total_num_models=0,
        structure_id="",
    )
    encoder.set_xtal_info(space_group="", unit_cell=None)

    encoder.set_header_info(
        r_free=None,
        r_work=None,
        resolution=None,
        title=None,
        deposition_date=None,
        release_date=None,
        experimental_methods=None,
    )

    node_ids = (
        df.model_id.astype(str)
        + ":"
        + df.chain_id
        + ":"
        + df.residue_name
        + ":"
        + df.residue_number.astype(str)
        + ":"
        + df.insertion.astype(str)
    )
    df["residue_id"] = node_ids
    # Tracks values to replace them at the end
    chains_per_model = []
    groups_per_chain = []

    # Iterate over each model
    for model_idx, model in enumerate(df.model_id.unique()):
        # Subset the DataFrame to the current model
        # Ang get the chains
        model_df = df.loc[df.model_id == model]
        chains = model_df.chain_id.unique()

        count_models += 1
        # Set the model info
        encoder.set_model_info(
            # model_id=model_idx, # According to mmtf-python this is meaningless
            model_id=model_idx,  # According to mmtf-python this is meaningless
            chain_count=0,  # Set to 0 here and changed later
        )
        # Iterate over chains in model
        for chain_id in chains:
            seqs = []
            seq = ""
            prev_res_type = ""
            prev_resname = ""
            first_chain = True

            # Subset the dataframe to the current chain
            chain_df = model_df[model_df.chain_id == chain_id]

            # Iterate over residues in chain
            residues = chain_df.residue_id.unique()
            for residue_id in residues:
                count_groups += 1
                # Subset the dataframe to the current residue
                residue_df = chain_df.loc[chain_df.residue_id == residue_id]
                # Get residue and entity types
                if all(residue_df.record_name == "ATOM"):
                    residue_type = "ATOM"
                    entity_type = "polymer"
                elif all(residue_df.residue_name == "HOH"):
                    residue_type = "HETATM"
                    entity_type = "water"
                else:
                    residue_type = "HETATM"
                    entity_type = "non-polymer"
                # Get the 3-letter residue name
                resname = residue_df.residue_name.iloc[0]

                # Check if the molecule changes within the chain
                # This will always increment for the first residue in a
                #  chain due to the starting values above
                # Checking for similar entities is non-trivial from the
                #  structure object so we treat each molecule as a separate
                #  entity
                if residue_type != prev_res_type or (
                    residue_type == "HETATM" and resname != prev_resname
                ):
                    encoder.set_entity_info(
                        chain_indices=[count_chains],
                        sequence="",  # Set to empty here and changed later
                        description="",
                        entity_type=entity_type,
                    )
                    encoder.set_chain_info(
                        chain_id=chain_id,
                        chain_name=(
                            "\x00" if len(chain_id.strip()) == 0 else chain_id
                        ),
                        num_groups=0,  # Set to 0 here and changed later
                    )
                    if count_chains > 0:
                        groups_per_chain.append(
                            count_groups - sum(groups_per_chain) - 1
                        )
                    if not first_chain:
                        seqs.append(seq)
                    first_chain = False
                    count_chains += 1
                    seq = ""

                if entity_type == "polymer":
                    seq += _seq1(
                        residue_df.residue_name.unique()[0],
                        charmap=protein_letters_3to1_extended,
                    )

                prev_res_type = residue_type
                prev_resname = resname

                group_type = (
                    "NON-POLYMER"
                    if residue_type == "HETATM"
                    else "L-PEPTIDE LINKING"
                )
                encoder.set_group_info(
                    group_name=residue_df.residue_name.unique()[0],
                    group_number=int(residue_df.residue_number.unique()[0]),
                    insertion_code=(
                        "\x00"
                        if residue_df.insertion.unique()[0] == ""
                        else residue_df.insertion.unique()[0]
                    ),
                    group_type=group_type,  # Hack to ensure we can re-parse.
                    atom_count=len(residue_df),
                    bond_count=0,
                    single_letter_code=_seq1(
                        df.residue_name.unique()[0],
                        charmap=protein_letters_3to1_extended,
                    ),
                    sequence_index=(
                        len(seq) - 1 if entity_type == "polymer" else -1
                    ),
                    secondary_structure_type=-1,
                )
                for row in residue_df.itertuples():
                    count_atoms += 1
                    encoder.set_atom_info(
                        atom_name=row.atom_name,
                        serial_number=row.atom_number,
                        alternative_location_id=(
                            "\x00" if row.alt_loc == "" else row.alt_loc
                        ),
                        x=row.x_coord,
                        y=row.y_coord,
                        z=row.z_coord,
                        occupancy=row.occupancy,
                        temperature_factor=row.b_factor,
                        element=row.element_symbol,
                        charge=0 if row.charge == "" else int(row.charge),
                    )
            seqs.append(seq)
            start_ind = len(encoder.entity_list) - len(seqs)
            for i, seq in enumerate(seqs):
                encoder.entity_list[start_ind + i]["sequence"] = seq

        chains_per_model.append(count_chains - sum(chains_per_model))
    groups_per_chain.append(count_groups - sum(groups_per_chain))

    encoder.chains_per_model = chains_per_model
    encoder.groups_per_chain = groups_per_chain
    encoder.num_atoms = count_atoms
    encoder.num_groups = count_groups
    encoder.num_chains = count_chains
    encoder.num_models = count_models
    encoder.finalize_structure()
    encoder.write_file(file_path)
