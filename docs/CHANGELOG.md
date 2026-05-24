# Release Notes ![](img/logos/3eiy_120.png)

- Supports `mol` files that have empty lines between blocks, (Via [Ruibin Liu](https://github.com/Ruibin-Liu) PR #[140](https://github.com/BioPandas/biopandas/pull/140#))

The CHANGELOG for the current development version is available at
[https://github.com/rasbt/biopandas/blob/main/docs/sources/CHANGELOG.md](https://github.com/rasbt/biopandas/blob/main/docs/sources/CHANGELOG.md).

### 0.6.1dev1 (UNRELEASED)
- Adds support for constructing a collection of PDB files and perform actions on them (Via [Julia K. Varga](https://github.com/gezmi)
- Adds the possibility of performing alignments, translation and rotation on the structure, via an extendable Align class (Via [Julia K. Varga](https://github.com/gezmi)
- TMAlign is installed alongside the package to perform the alignments (Via [Julia K. Varga](https://github.com/gezmi)
- Added feature to write out mmcif files
- Fixed erronoues links in documentation, together with automatic generation of said documents (Via [Julia K. Varga](
- Dev: switched testing framework entirely to pytest. Drops nose dependency due to version conflicts with Python 3.12 (`nose`) and 3.8 (`nose`)

### 0.5.1 (01/08/2024)

- Fix: improves support for writing PDBs with `OTHERS` records to stream. PR [#149](https://github.com/BioPandas/biopandas/pull/149). Addresses issue [#141](https://github.com/BioPandas/biopandas/issues/141).
- Feature: added method to `PandasMmcif` that allow to select by model ids. PR #[145](https://github.com/BioPandas/biopandas/pull/145)
- Dev: switched testing framework entirely to pytest. Drops nose dependency due to version conflicts with Python 3.12 (`nose`) and 3.8 (`nose`) PR #[146](https://github.com/BioPandas/biopandas/pull/146)
- Dev: adds GitHub actions-based CI. PR [#149](https://github.com/BioPandas/biopandas/pull/149).
- Avoid inclusion of test scripts and test data in the PyPI release of the Biopandas library. PR #[148](https://github.com/BioPandas/biopandas/pull/148). Addresses issue [#147](https://github.com/BioPandas/biopandas/issues/147)

### 0.5.0dev1 (31/7/2023)

- Implement add_remark for PandasPdb, (Via [Anton Bushuiev](https://github.com/anton-bushuiev) PR #[129](https://github.com/BioPandas/biopandas/pull/129))
- B_factor shifting one white space issue fix. (Via [Zehra Sarica](https://github.com/zehraacarsarica), PR #[134](https://github.com/BioPandas/biopandas/pull/134))
- Adds support for pathlib. (Via [Anton Bushuiev](https://github.com/anton-bushuiev), PR #[128](https://github.com/BioPandas/biopandas/pull/128))
- Adds support for reading Gzipped MMTF files. (Via [Arian Jamasb](https://github.com/a-r-j), PR #[123](https://github.com/rasbt/biopandas/pull/123/files))
- Improves reliability of parsing polymer/non-polymer entities in MMTF parsing. (Via [Arian Jamasb](https://github.com/a-r-j), PR #[123](https://github.com/rasbt/biopandas/pull/123/files))
- Improves reliability of parsing multicharacter chain IDs from MMTF files. (Via [Arian Jamasb](https://github.com/a-r-j), PR #[123](https://github.com/rasbt/biopandas/pull/123/files))
- Replaces null terminator chars in parsed MMTF dataframe with the empty string. (Via [Arian Jamasb](https://github.com/a-r-j), PR #[123](https://github.com/rasbt/biopandas/pull/123/files))

### 0.5.0dev0 (3/4/2023)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.5.0.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.5.0.tar.gz)

##### New Features

- Added a new `PandasMmcif.convert_to_pandas_pdb()` class that converts the mmCIF file into a PDB structure. (Via [Arian Jamasb](https://github.com/a-r-j), PR #[107](https://github.com/rasbt/biopandas/pull/107/files))
- Added ability to write PDBs to a filestream (Via [Arian Jamasb](https://github.com/a-r-j), PR #[107](https://github.com/rasbt/biopandas/pull/107/files))
- Adds a method `PandasPdb.gyradius` to calculate the radius of gyration of a molecule (via [goniochromatic](https://github.com/github.com/goniochromatic/), PR #[117](https://github.com/rasbt/biopandas/pull/117/files))
- Adds MMTF export & improves MMTF parsing robustness (via [Arian Jamasb](https://github.com/a-r-j), PR #[119](https://github.com/rasbt/biopandas/pull/119/files))
- Adds support for parsing [MMTF](https://mmtf.rcsb.org/) files. (via [Arian Jamasb](https://github.com/a-r-j), PR #[111](https://github.com/rasbt/biopandas/pull/111/files))
- Adds primitive functions for parsing PDB, mmCIF, and MMTF into dataframes. (via [Arian Jamasb](https://github.com/a-r-j), PR #[111](https://github.com/rasbt/biopandas/pull/111/files))
- Added support for [AlphaFolds 200M+ structures](https://www.deepmind.com/blog/alphafold-reveals-the-structure-of-the-protein-universe) via `PandasMmcif().fetch_mmcif(uniprot_id='Q5VSL9', source='alphafold2-v3')` and `PandasPdb().fetch_pdb(uniprot_id='Q5VSL9', source='alphafold2-v3')`. (Via [Arian Jamasb](https://github.com/a-r-j), PR #[102](https://github.com/rasbt/biopandas/pull/102/files))

##### Bug Fixes

- Fix the `return` statement in `PandasPdb.to_pdb_stream()` to return `output` instead of `output.seek(0)`. (via [goniochromatic](https://github.com/github.com/goniochromatic/), PR #[116](https://github.com/rasbt/biopandas/pull/116/files))
- Change the `records` default argument in `PandasPdb.to_pdb_stream()` to be immutable. (via [goniochromatic](https://github.com/github.com/goniochromatic/), PR #[116](https://github.com/rasbt/biopandas/pull/116/files))
- Fix some typos and general style issues. (via [goniochromatic](https://github.com/github.com/goniochromatic/), PR #[116](https://github.com/rasbt/biopandas/pull/116/files))
- Fix link for "How to contribute" in `README.md`. (via [goniochromatic](https://github.com/github.com/goniochromatic/), PR #[116](https://github.com/rasbt/biopandas/pull/116/files))

### 0.4.1 (05-13-2022)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.4.1.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.4.1.tar.gz)

##### Changes

- Remove walrus operator for Python 3.7 compatibility.

### 0.4.0 (05-11-2022)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.4.0.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.4.0.tar.gz)

##### New Features

- Adds support for extracting structures from PDB files containing multiple models. See the [documentation](http://rasbt.github.io/biopandas/tutorials/Working_with_PDB_Structures_in_DataFrames/#working-with-pdbs-containing-multiple-models
) for details.  (via [Arian Jamasb](https://github.com/a-r-j), PR #[101](https://github.com/rasbt/biopandas/pull/101/files)).

- Adds support for fetching mmCIF (`PandasMmcif().fetch_mmcif(uniprot_id='Q5VSL9', source='alphafold2-v2')`) and PDB structures (e.g., `PandasPdb().fetch_pdb(uniprot_id='Q5VSL9', source="alphafold2-v2")`)  (via [Arian Jamasb](https://github.com/a-r-j), PR #[102](https://github.com/rasbt/biopandas/pull/102/files)).

##### Changes

- Instead of raising a warning when no ATOM entries are loaded, raise the warning only when neither ATOM nor HETAM entries are loaded.

##### Bug Fixes

- None

### 0.3.0 (04-06-2022)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.3.0.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.3.0.tar.gz)

##### New Features

- Adds support for parsing mmCIF protein structure files (via [Arian Jamasb](https://github.com/a-r-j), PR #[94](https://github.com/rasbt/biopandas/pull/94/files))

##### Changes

- -

##### Bug Fixes

- Fixes a bug where coordinates with more than 4 digits before the decimal point caused a column shift when saving a PDB file. (via  PR #[90](https://github.com/rasbt/biopandas/pull/90/files))
- Fixes a bug where the invert parameter in get_carbon was selecting the wrong case. (via [Arian Jamasb](https://github.com/a-r-j) PR #[96](https://github.com/rasbt/biopandas/pull/96/files))

### 0.2.9 (08-30-2021)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.9.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.9.tar.gz)

##### New Features

- -

##### Changes

- Now also allow `.ent` and `.ent.gz` file endings for PDB files.  (via PR #[82](https://github.com/rasbt/biopandas/pull/82/files)
- Added Python 3.8 and 3.9 to setup.py in order to support these versions via conda-forge. (via PR #[87](https://github.com/rasbt/biopandas/pull/87/files)

##### Bug Fixes

- -

### 0.2.8 (03-30-2021)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.8.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.8.tar.gz)

##### New Features

- A `PandasPdb.read_pdb_from_list` method was added analogous to the existing `PandasMol2.read_mol2_from_list` (via PR #[72](https://github.com/rasbt/biopandas/pull/72/files) by [dominiquesydow](https://github.com/dominiquesydow))

##### Changes

- `ValueError` raising and improved file format error messages for `read_pdb` and `read_mol2` functionality. (via PR #[73](https://github.com/rasbt/biopandas/pull/73/files) by [dominiquesydow](https://github.com/dominiquesydow))

##### Bug Fixes

- -

### 0.2.7 (08-04-2020)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.7.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.7.tar.gz)

##### New Features

- -

##### Changes

- -

##### Bug Fixes

- Fix Manifest file to include license file in the PyPI tar.gz file so that BioPandas can be packaged by conda-forge.

### 0.2.6 (08-03-2020)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.6.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.6.tar.gz)

##### New Features

- -

##### Changes

- Uses more modern `https` queries for the RCSB server via the `fetch_pdb` function.
- Updates the documentation (incl. a code of conduct)

##### Bug Fixes

- -

### 0.2.5 (07-09-2019)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.5.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.5.tar.gz)

##### New Features

- -

##### Changes

- -

##### Bug Fixes

- The `PandasPdb.amino3to1` method now also considers insertion codes when converting the amino acid codes; before, inserted amino acides were skipped.

### 0.2.4 (02-05-2019)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.4.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.4.tar.gz)

##### New Features

- -

##### Changes

- Minor adjustments to support to address deprecation warnings in pandas >= 23.0

##### Bug Fixes

- -

### 0.2.3 (03-29-2018)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.3.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.3.tar.gz)

##### New Features

- -

##### Changes

- `PandasMol2.distance_df` was added as a static method that allows distance computations based for external data frames with its behavior otherwise similar to `PandasMol2.distance`.
- `PandasPdb.distance_df` was added as a static method that allows distance computations based for external data frames with its behavior otherwise similar to `PandasPdb.distance`.
- `PandasPdb.distance` now supports multiple record sections to be considered (e.g., `records=('ATOM', 'HETATM')` to include both protein and ligand in a query. Now also defaults to `records=('ATOM', 'HETATM')` for concistency with the impute method.
- `PandasPdb.get(...)` now supports external data frames and lets the user specify the record section to be considered (e.g., `records=('ATOM', 'HETATM')` to include both protein and ligand in a query. Now also defaults to `records=('ATOM', 'HETATM')` for concistency with the impute method.
- The `section` parameter of `PandasPdb.impute_element(...)` was renamed to `records` for API consistency.

##### Bug Fixes

-

### 0.2.2 (06-07-2017)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.2.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.2.tar.gz)

##### New Features

- -

##### Changes

- Raises a meaningful error message if attempting to overwrite the `df` attributes of `PandasMol2` and `PandasPdb` directly.
- Added `PandasPdb.pdb_path` and `PandasMol2.mol2_path` attributes that store the location of the data file last read.

##### Bug Fixes

- The `rmsd` methods of `PandasMol2` and `PandasPdb` don't return a NaN anymore if the array indices of to structures are different.

### 0.2.1  (2017-05-11)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.1.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.1.tar.gz)

##### New Features

- -

##### Changes

- The `amino3to1` method of `biopandas.pdb.PandasPDB` objects now returns a pandas `DataFrame` instead of a pandas `Series` object. The returned data frame has two columns, `'chain_id'` and `'residue_name'`, where the former contains the chain ID of the amino acid and the latter contains the 1-letter amino acid code, respectively.
- Significant speed improvements of the `distance` method of both `PandasPdb` and `PandasMol2` (now about 300 percent faster than previously).

##### Bug Fixes

- The `amino3to1` method of `biopandas.pdb.PandasPDB` objects now handles multi-chain proteins correctly.
- The `amino3to1` method of `biopandas.pdb.PandasPDB` objects now also works as expected if the `'ATOM'` entry DataFrame contains disordered DataFrame indices or duplicate DataFrame index values.

### 0.2.0 (2017-04-02)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/archive/v0.2.0.zip)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/archive/v0.2.0.tar.gz)

##### New Features

- Added an `amino3to1` method to `PandasPdb` data frames to convert 3-amino acid letter codes to 1-letter codes.
- Added a `distance` method to `PandasPdb` data frames to compute the Euclidean distance between atoms and a reference point.
- Added the `PandasMol2` class for working with Tripos MOL2 files in pandas DataFrames.

##### Changes

- `PandasPDB` was renamed to `PandasPdb`.
- Raises a warning if `PandasPdb` is written to PDB and ATOM and HETAM section contains unexpected columns; these columns will now be skipped.

##### Bug Fixes

- -

### 0.1.5 (2016-11-19)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/releases/tag/v0.1.5)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/releases/tag/v0.1.5.tar.gz)

##### New Features

- Added an `impute_element` method to `PandasPDB` objects to infer the Element Symbol from the Atom Name column.
- Added two new selection types for `PandasPDB` ATOM and HETATM coordinate sections: `'heavy'` and `'carbon'`.

##### Changes

- Include test data in the PyPI package; add install_requires for pandas.
- The `'hydrogen'` atom selection in `PandasPDB` methods is now based on the element type instead of the atom name.
- By default, the RMSD is now computed on all atoms unless a specific selection is defined.

##### Bug Fixes

- -

### 0.1.4 (2015-11-24)

##### Downloads

- [Source code (zip)](https://github.com/rasbt/biopandas/releases/tag/v0.1.4)
- [Source code (tar.gz)](https://github.com/rasbt/biopandas/releases/tag/v0.1.4.tar.gz)

##### New Features

-

##### Changes

- Needed to bump the version number due to a bug in the PyPI setup.py script.
- Support for the old pandas sorting syntax (`DataFrame.sort` vs `DataFrame.sort_values`) incl. DeprecationWarning.

##### Bug Fixes

-

### 0.1.3 (2015-11-23)

##### New Features

-

##### Changes

-

##### Bug Fixes

- Exception handling in tests if PDB goes down (which just happened).
- Added a separate ANISOU engine to handle those records correctly.

### 0.1.2 (2015-11-23)

- First Release.
