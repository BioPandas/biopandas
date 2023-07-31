# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

import pandas as pd

amino3to1dict = {
    "ASH": "A",
    "ALA": "A",
    "CYX": "C",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "HID": "H",
    "HIE": "H",
    "HIP": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "MSE": "M",
    "ASN": "N",
    "PYL": "O",
    "HYP": "P",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "SEL": "U",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}

pdb_df_columns = {
    "record_name",
    "atom_number",
    "blank_1",
    "atom_name",
    "alt_loc",
    "residue_name",
    "blank_2",
    "chain_id",
    "residue_number",
    "insertion",
    "blank_3",
    "x_coord",
    "y_coord",
    "z_coord",
    "occupancy",
    "b_factor",
    "blank_4",
    "segment_id",
    "element_symbol",
    "charge",
}

pdb_atomdict = [
    {"id": "record_name", "line": [0, 6], "type": str, "strf": lambda x: "%-6s" % x},
    {
        "id": "atom_number",
        "line": [6, 11],
        "type": int,
        "strf": lambda x: "%+5s" % str(x),
    },
    {"id": "blank_1", "line": [11, 12], "type": str, "strf": lambda x: "%-1s" % x},
    {
        "id": "atom_name",
        "line": [12, 16],
        "type": str,
        "strf": lambda x: " %-3s" % x if len(x) < 4 else "%-4s" % x,
    },
    {"id": "alt_loc", "line": [16, 17], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "residue_name", "line": [17, 20], "type": str, "strf": lambda x: "%+3s" % x},
    {"id": "blank_2", "line": [20, 21], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "chain_id", "line": [21, 22], "type": str, "strf": lambda x: "%-1s" % x},
    {
        "id": "residue_number",
        "line": [22, 26],
        "type": int,
        "strf": lambda x: "%+4s" % str(x),
    },
    {"id": "insertion", "line": [26, 27], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "blank_3", "line": [27, 30], "type": str, "strf": lambda x: "%-3s" % x},
    {
        "id": "x_coord",
        "line": [30, 38],
        "type": float,
        "strf": lambda x: ("%+8.3f" % x).replace("+", " "),
    },
    {
        "id": "y_coord",
        "line": [38, 46],
        "type": float,
        "strf": lambda x: ("%+8.3f" % x).replace("+", " "),
    },
    {
        "id": "z_coord",
        "line": [46, 54],
        "type": float,
        "strf": lambda x: ("%+8.3f" % x).replace("+", " "),
    },
    {
        "id": "occupancy",
        "line": [54, 60],
        "type": float,
        "strf": lambda x: ("%+6.2f" % x).replace("+", " "),
    },
    {
        "id": "b_factor",
        "line": [60, 66],
        "type": float,
        "strf": lambda x: ("%+6.2f" % x).replace("+", " ") if len(str(int(x))) < 3 else ("%+6.2f" % x).replace("+", ""),
    },
    {"id": "blank_4", "line": [66, 72], "type": str, "strf": lambda x: "%-7s" % x},
    {"id": "segment_id", "line": [72, 76], "type": str, "strf": lambda x: "%-3s" % x},
    {
        "id": "element_symbol",
        "line": [76, 78],
        "type": str,
        "strf": lambda x: "%+2s" % x,
    },
    {
        "id": "charge",
        "line": [78, 80],
        "type": float,
        "strf": lambda x: (("%+2.1f" % x).replace("+", " ") if pd.notnull(x) else ""),
    },
]


pdb_anisoudict = [
    {"id": "record_name", "line": [0, 6], "type": str, "strf": lambda x: "%-6s" % x},
    {
        "id": "atom_number",
        "line": [6, 11],
        "type": int,
        "strf": lambda x: "%+5s" % str(x),
    },
    {"id": "blank_1", "line": [11, 12], "type": str, "strf": lambda x: "%-1s" % x},
    {
        "id": "atom_name",
        "line": [12, 16],
        "type": str,
        "strf": lambda x: (" %-3s" % x if len(x) < 4 else "%-4s" % x),
    },
    {"id": "alt_loc", "line": [16, 17], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "residue_name", "line": [17, 20], "type": str, "strf": lambda x: "%+3s" % x},
    {"id": "blank_2", "line": [20, 21], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "chain_id", "line": [21, 22], "type": str, "strf": lambda x: "%-1s" % x},
    {
        "id": "residue_number",
        "line": [22, 26],
        "type": int,
        "strf": lambda x: "%+4s" % str(x),
    },
    {"id": "insertion", "line": [26, 27], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "blank_3", "line": [27, 28], "type": str, "strf": lambda x: "%-1s" % x},
    {"id": "U(1,1)", "line": [28, 35], "type": int, "strf": lambda x: "%+7s" % str(x)},
    {"id": "U(2,2)", "line": [35, 42], "type": int, "strf": lambda x: "%+7s" % str(x)},
    {"id": "U(3,3)", "line": [42, 49], "type": int, "strf": lambda x: "%+7s" % str(x)},
    {"id": "U(1,2)", "line": [49, 56], "type": int, "strf": lambda x: "%+7s" % str(x)},
    {"id": "U(1,3)", "line": [56, 63], "type": int, "strf": lambda x: "%+7s" % str(x)},
    {"id": "U(2,3)", "line": [63, 70], "type": int, "strf": lambda x: "%+7s" % str(x)},
    {"id": "blank_4", "line": [70, 76], "type": str, "strf": lambda x: "%+6s" % x},
    {
        "id": "element_symbol",
        "line": [76, 78],
        "type": str,
        "strf": lambda x: "%+2s" % x,
    },
    {
        "id": "charge",
        "line": [78, 80],
        "type": float,
        "strf": lambda x: (("%+2.1f" % x).replace("+", " ") if pd.notnull(x) else ""),
    },
]

pdb_otherdict = [
    {
        "id": "record_name",
        "line": [0, 6],
        "type": str,
        "strf": lambda x: "%s%s" % (x, " " * (6 - len(x))),
    },
    {"id": "entry", "line": [6, -2], "type": str, "strf": lambda x: x.rstrip()},
]

pdb_records = {
    "ATOM": pdb_atomdict,
    "HETATM": pdb_atomdict,
    "ANISOU": pdb_anisoudict,
    "OTHERS": pdb_otherdict,
}
