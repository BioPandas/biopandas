"""Utility function for reading Tripos MOL2 files from files"""

# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import sys

if sys.version_info >= (3, 9):
    import importlib.resources as pkg_resources
else:
    import importlib_resources as pkg_resources

import tests.mol2.data
from biopandas.mol2 import PandasMol2
from biopandas.mol2.mol2_io import split_multimol2
from tests.testutils import assert_raises

TEST_DATA = pkg_resources.files(tests.mol2.data)


def test_read_mol2():
    data_path_1 = str(TEST_DATA.joinpath("40_mol2_files.mol2"))
    data_path_2 = str(TEST_DATA.joinpath("40_mol2_files.mol2.gz"))
    data_path_3 = str(TEST_DATA.joinpath("empty_line.mol2"))

    for data_path in (data_path_1, data_path_2):
        pdmol = PandasMol2().read_mol2(data_path)
        assert pdmol.df.shape == (65, 9)
        assert pdmol.code == "ZINC38611810"

        expect = [
            "atom_id",
            "atom_name",
            "x",
            "y",
            "z",
            "atom_type",
            "subst_id",
            "subst_name",
            "charge",
        ]
        assert expect == list(pdmol.df.columns)
        assert len(pdmol.mol2_text) == 6469
        assert pdmol.mol2_path == data_path

    pdmol = PandasMol2().read_mol2(data_path_3)
    assert pdmol.df.shape == (59, 9)


def test_read_mol2_from_list():
    data_path = str(TEST_DATA.joinpath("40_mol2_files.mol2"))
    mol2 = next(split_multimol2(data_path))

    pdmol = PandasMol2().read_mol2_from_list(
        mol2_lines=mol2[1], mol2_code=mol2[0]
    )
    assert pdmol.df.shape == (65, 9)
    assert pdmol.code == "ZINC38611810"


def test_rmsd():
    data_path_1 = str(TEST_DATA.joinpath("1b5e_1.mol2"))
    data_path_2 = str(TEST_DATA.joinpath("1b5e_2.mol2"))

    pdmol_1 = PandasMol2().read_mol2(data_path_1)
    pdmol_2 = PandasMol2().read_mol2(data_path_2)

    assert pdmol_1.rmsd(pdmol_1.df, pdmol_2.df, heavy_only=False) == 1.5523
    # Divides by the number of heavy atoms (20), not by the total atom
    # count (32) that includes the hydrogens filtered out above.
    assert pdmol_1.rmsd(pdmol_1.df, pdmol_2.df) == 1.4684


def test_rmsd_heavy_only_normalizes_by_compared_atoms():
    """`heavy_only=True` must normalize by the atoms it actually compared.

    The sum of squared deviations is accumulated over the hydrogen-filtered
    frames, so dividing it by the unfiltered atom count understates the
    result by a factor of sqrt(n_total / n_heavy).
    """
    data_path_1 = str(TEST_DATA.joinpath("1b5e_1.mol2"))
    data_path_2 = str(TEST_DATA.joinpath("1b5e_2.mol2"))

    df1 = PandasMol2().read_mol2(data_path_1).df
    df2 = PandasMol2().read_mol2(data_path_2).df

    # The fixtures must contain hydrogens, otherwise this test cannot detect
    # the mis-normalization.
    n_total = df1.shape[0]
    heavy1 = df1[df1["atom_type"] != "H"]
    heavy2 = df2[df2["atom_type"] != "H"]
    n_heavy = heavy1.shape[0]
    assert n_heavy < n_total

    sq = (
        (heavy1["x"].values - heavy2["x"].values) ** 2
        + (heavy1["y"].values - heavy2["y"].values) ** 2
        + (heavy1["z"].values - heavy2["z"].values) ** 2
    )
    expected = round((sq.sum() / n_heavy) ** 0.5, 4)

    actual = PandasMol2.rmsd(df1, df2)
    assert actual == expected
    # Guard against the value that the unfiltered divisor would produce.
    assert actual != round((sq.sum() / n_total) ** 0.5, 4)


def test_distance():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))

    pdmol = PandasMol2().read_mol2(data_path)
    assert round(pdmol.distance().values[0], 3) == 31.185


def test_distance_external_df():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))

    pdmol = PandasMol2().read_mol2(data_path)
    new_df = pdmol.df.iloc[1:, :].copy()
    assert round(PandasMol2.distance_df(df=new_df).values[0], 3) == 31.165


def test_overwrite_df():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))
    pdmol = PandasMol2().read_mol2(data_path)

    def overwrite():
        pdmol.df = pdmol.df[(pdmol.df["atom_type"] != "H")]

    expect = (
        "Please use `PandasMol2._df = ... `"
        " instead\nof `PandasMol2.df = ... `"
        " if you are sure that\nyou want"
        " to overwrite the `df` attribute."
    )

    assert_raises(AttributeError, expect, overwrite)


def test__get_atomsection_raises():
    """Test if ValueError is raised if input list is not in the mol2 format."""

    expect = (
        "Structural data could not be loaded. "
        "Is the input file/text in the mol2 format?"
    )

    def run_code():
        PandasMol2()._get_atomsection(["", ""])

    assert_raises(ValueError, expect, run_code)


def test_parse_bonds():
    data_path = str(TEST_DATA.joinpath("1b5e_1.mol2"))
    pdmol = PandasMol2().read_mol2(data_path)

    # Bonds are only parsed on request.
    assert pdmol.df_bonds is None

    bonds = pdmol.parse_bonds()
    assert bonds is pdmol.df_bonds
    assert bonds.shape == (33, 5)
    assert list(bonds.columns) == [
        "bond_id",
        "origin_atom_id",
        "target_atom_id",
        "bond_type",
        "status_bits",
    ]
    for col in ("bond_id", "origin_atom_id", "target_atom_id"):
        assert bonds[col].dtype == "int64"
    assert bonds.iloc[0].tolist() == [1, 1, 2, "2", ""]
    # Non-numeric Tripos bond types are kept as strings.
    assert bonds["bond_type"].value_counts().to_dict() == {
        "1": 28,
        "2": 4,
        "am": 1,
    }
    # Every bond connects atoms of the ATOM section.
    assert set(bonds["origin_atom_id"]) <= set(pdmol.df["atom_id"])
    assert set(bonds["target_atom_id"]) <= set(pdmol.df["atom_id"])


def test_parse_bonds_matches_header_counts():
    """Bond counts match the MOLECULE header for every test molecule,
    including gzipped files and files with CRLF line endings."""
    for fname in (
        "40_mol2_files.mol2",
        "40_mol2_files.mol2.gz",
        "1b5e_2.mol2",
        "empty_line.mol2",
    ):
        data_path = str(TEST_DATA.joinpath(fname))
        for code, lines in split_multimol2(data_path):
            pdmol = PandasMol2().read_mol2_from_list(lines, code)
            header = [
                s for s in pdmol.mol2_text.splitlines() if s.strip()
            ]
            start = header.index(
                next(s for s in header if s.startswith("@<TRIPOS>MOLECULE"))
            )
            num_bonds = int(header[start + 2].split()[1])
            assert len(pdmol.parse_bonds()) == num_bonds, (fname, code)


def test_parse_bonds_status_bits_and_end_of_file():
    """Status bits are parsed, blank and comment lines are skipped, and a
    BOND section at the end of the file is read completely."""
    mol2_lines = [
        "@<TRIPOS>MOLECULE\n",
        "ethene\n",
        "2 1 0 0 0\n",
        "SMALL\n",
        "NO_CHARGES\n",
        "\n",
        "@<TRIPOS>ATOM\n",
        "1 C1 0.0000 0.0000 0.0000 C.2 1 ETH 0.0000\n",
        "2 C2 1.3300 0.0000 0.0000 C.2 1 ETH 0.0000\n",
        "\n",
        "@<TRIPOS>BOND\n",
        "# a comment\n",
        "\n",
        "     1     1     2 2 BACKBONE|DICT",
    ]
    pdmol = PandasMol2().read_mol2_from_list(mol2_lines, "ethene")
    bonds = pdmol.parse_bonds()
    assert bonds.shape == (1, 5)
    assert bonds.iloc[0].tolist() == [1, 1, 2, "2", "BACKBONE|DICT"]


def test_parse_bonds_without_bond_section():
    mol2_lines = [
        "@<TRIPOS>MOLECULE\n",
        "sodium\n",
        "1 0 0 0 0\n",
        "SMALL\n",
        "USER_CHARGES\n",
        "\n",
        "@<TRIPOS>ATOM\n",
        "1 Na 0.0000 0.0000 0.0000 Na 1 NA 1.0000\n",
        "\n",
    ]
    pdmol = PandasMol2().read_mol2_from_list(mol2_lines, "sodium")
    bonds = pdmol.parse_bonds()
    assert bonds.shape == (0, 5)
    assert bonds["bond_id"].dtype == "int64"


def test_parse_bonds_resets_on_new_read():
    pdmol = PandasMol2().read_mol2(str(TEST_DATA.joinpath("1b5e_1.mol2")))
    pdmol.parse_bonds()
    pdmol.read_mol2(str(TEST_DATA.joinpath("40_mol2_files.mol2")))
    assert pdmol.df_bonds is None


def test_parse_bonds_raises():
    expect = (
        "No Mol2 structure loaded. Please read a Mol2 file "
        "via `read_mol2` or `read_mol2_from_list` first."
    )
    assert_raises(AttributeError, expect, PandasMol2().parse_bonds)

    mol2_lines = [
        "@<TRIPOS>MOLECULE\n",
        "broken\n",
        "@<TRIPOS>ATOM\n",
        "1 C1 0.0000 0.0000 0.0000 C.3 1 UNK 0.0000\n",
        "\n",
        "@<TRIPOS>BOND\n",
        "1 1 2\n",
    ]
    pdmol = PandasMol2().read_mol2_from_list(mol2_lines, "broken")
    expect = (
        "Could not parse BOND record '1 1 2'. Expected at least 4 "
        "fields: bond_id origin_atom_id target_atom_id "
        "bond_type [status_bits]"
    )
    assert_raises(ValueError, expect, pdmol.parse_bonds)
