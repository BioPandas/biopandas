"""Tests for observed protein sequence export."""

from copy import deepcopy
from pathlib import Path

import pandas as pd
import pytest

from biopandas.pdb import PandasPdb


DATA = Path(__file__).parent / "data"


def protein(models=None):
    if models is None:
        models = [[("B", "MET"), ("B", "LYS"), ("A", "GLY")]]
    lines = []
    for model, residues in enumerate(models, 1):
        if len(models) > 1:
            lines.append(f"MODEL     {model:4d}\n")
        for serial, (chain, residue) in enumerate(residues, 1):
            lines.append(
                f"ATOM  {serial:5d}  CA  {residue:3s} "
                f"{chain:1s}{serial:4d}    "
                f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}"
                "           C  \n"
            )
        if len(models) > 1:
            lines.append("ENDMDL\n")
    return PandasPdb().read_pdb_from_list(lines)


def test_chain_order_and_exact_output():
    assert protein().to_fasta(model_index=1) == (
        ">model_1_chain_B\nMK\n>model_1_chain_A\nG\n"
    )


@pytest.mark.parametrize("width,sequence", [(1, "M\nK"), (2, "MK"), (3, "MK")])
def test_wrapping(width, sequence):
    pdb = protein([[("A", "MET"), ("A", "LYS")]])
    assert pdb.to_fasta(model_index=1, line_width=width) == (
        f">model_1_chain_A\n{sequence}\n"
    )


@pytest.mark.parametrize("chain", ["", " "])
def test_blank_chain(chain):
    pdb = protein([[("A", "MET")]])
    pdb.df["ATOM"]["chain_id"] = chain
    assert pdb.to_fasta(model_index=1) == ">model_1_chain_blank\nM\n"


@pytest.mark.parametrize("chain", ["A B", "\t", "\n", "A\x00", "  ", None])
def test_invalid_chain(chain):
    pdb = protein()
    pdb.df["ATOM"]["chain_id"] = chain
    with pytest.raises(ValueError, match="Chain identifiers"):
        pdb.to_fasta(model_index=1)


def test_translation():
    pdb = protein([[("A", "UNK"), ("A", "MSE"), ("A", "HYP")]])
    assert pdb.to_fasta(model_index=1) == ">model_1_chain_A\nXMP\n"
    assert pdb.to_fasta(model_index=1, fillna="Z") == ">model_1_chain_A\nZMP\n"


def test_insertion_codes():
    pdb = PandasPdb().read_pdb(DATA / "2d7t.pdb")
    sequence = "".join(pdb.to_fasta(model_index=1).splitlines()[1:])
    assert sequence[50:60] == "INPKSGDTNY"


def test_model_selection():
    pdb = protein([[("A", "MET")], [("A", "GLY")]])
    assert pdb.to_fasta(model_index=1) == ">model_1_chain_A\nM\n"
    assert pdb.to_fasta(model_index=2) == ">model_2_chain_A\nG\n"


def test_absent_model():
    with pytest.raises(ValueError, match="absent"):
        protein().to_fasta(model_index=2)


@pytest.mark.parametrize("name,value", [
    ("model_index", True), ("model_index", "1"), ("model_index", 1.0),
    ("line_width", False), ("line_width", "80"), ("line_width", 1.5),
])
def test_integer_arguments(name, value):
    arguments = {"model_index": 1, name: value}
    with pytest.raises(TypeError, match="integer"):
        protein().to_fasta(**arguments)


@pytest.mark.parametrize("width", [0, -1])
def test_invalid_width(width):
    with pytest.raises(ValueError, match="positive"):
        protein().to_fasta(model_index=1, line_width=width)


@pytest.mark.parametrize("fillna", [None, "", "XX", "x", "?", "É"])
def test_invalid_fillna(fillna):
    with pytest.raises(ValueError, match="uppercase ASCII"):
        protein().to_fasta(model_index=1, fillna=fillna)


def test_unloaded_and_required_model():
    with pytest.raises(ValueError, match="No PDB structure"):
        PandasPdb().to_fasta(model_index=1)
    with pytest.raises(TypeError):
        protein().to_fasta()


def test_empty_selection(tmp_path):
    pdb = protein()
    pdb.df["ATOM"] = pdb.df["ATOM"].iloc[:0]
    path = tmp_path / "empty.fasta"
    assert pdb.to_fasta(path, model_index=1) == ""
    assert path.read_bytes() == b""


def test_hetatm_excluded():
    pdb = protein()
    pdb.df["HETATM"] = pdb.df["ATOM"].copy()
    pdb.df["HETATM"]["residue_name"] = "ALA"
    assert pdb.to_fasta(model_index=1) == (
        ">model_1_chain_B\nMK\n>model_1_chain_A\nG\n"
    )


@pytest.mark.parametrize("as_string", [False, True])
def test_file_output(tmp_path, as_string):
    path = tmp_path / "chains.fasta"
    path.write_text("old content")
    fasta = protein().to_fasta(str(path) if as_string else path, model_index=1)
    assert path.read_bytes() == fasta.encode("utf-8")
    assert fasta.endswith("\n") and not fasta.endswith("\n\n")


def test_source_unchanged():
    pdb = protein([[("A", "MET")], [("A", "GLY")]])
    before = deepcopy(pdb.df)
    text = pdb.pdb_text
    pdb.to_fasta(model_index=2)
    assert pdb.df.keys() == before.keys()
    for record in before:
        pd.testing.assert_frame_equal(before[record], pdb.df[record])
    assert pdb.pdb_text == text


@pytest.mark.parametrize("invalid_chain", [False, True])
def test_validation_preserves_destination(tmp_path, invalid_chain):
    path = tmp_path / "existing.fasta"
    path.write_text("keep this")
    pdb = protein()
    if invalid_chain:
        pdb.df["ATOM"]["chain_id"] = "A\n"
    with pytest.raises(ValueError):
        pdb.to_fasta(path, model_index=1 if invalid_chain else 99)
    assert path.read_text() == "keep this"
