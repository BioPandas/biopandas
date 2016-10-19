# Changelog ![](img/logos/3eiy_120.png)

#### 0.1.5

- Include test data in the PyPI package; add install_requires for pandas
- Added an `impute_element` method to `PandasPDB` objects to infer the Element Symbol from the Atom Name column
- Added two new selection types for `PandasPDB` ATOM and HETATM coordinate sections: `'heavy'` and `'carbon'`
- The `'hydrogen'` atom selection in `PandasPDB` methods is now based on the element type instead of the atom name
- By default, the RMSD is now computed on all atoms unless a specific selection is defined

#### 0.1.4 (2015-11-24)

- Needed to bump the version number due to a bug in the PyPI setup.py script
- Support for the old pandas sorting syntax (`DataFrame.sort` vs `DataFrame.sort_values`) incl. DeprecationWarning

#### 0.1.3 (2015-11-23)

- Exception handling in tests if PDB goes down (which just happened)
- Added a separate ANISOU engine to handle those records correctly


#### 0.1.2 (2015-11-23)

- First Release
