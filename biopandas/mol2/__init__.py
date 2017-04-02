# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas

"""
BioPandas module for working with TRIPOS MOL2
files in pandas DataFrames.
"""

from .pandas_mol2 import PandasMol2
from .mol2_io import split_multimol2

__all__ = ["PandasMol2", "split_multimol2"]
