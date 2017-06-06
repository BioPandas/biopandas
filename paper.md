---
title: 'BioPandas: Working with molecular structures in pandas DataFrames'
tags:
  - bioinformatics
  - computational biology
  - protein structure analysis
  - protein-ligand docking
  - virtual screening
authors:
 - name: Sebastian Raschka
   orcid: 0000-0001-6989-4493
   affiliation: 1
affiliations:
 - name: Michigan State University, East-Lansing, USA
   index: 1
date: 31 May 2017
---



# Summary

BioPandas is a Python library that reads molecular structures from 3D-coordinate files, such as PDB (Bernstein *and others* 1977) and MOL2, into pandas DataFrames (McKinney 2010) for convenient data analysis and data mining related tasks.

In addition to parsing protein and small molecule data into a data frame format, BioPandas provides additional utility functions for structure analysis. These functions include common computations such as computing the root-mean-squared-deviation between structures and converting protein structures into primary amino acid sequence formats. 

Furthermore, useful small-molecule related functions are provided for reading and parsing millions of small molecule structures (from multi-MOL2 files) fast and efficiently in virtual screening applications. Inbuilt functions for filtering molecules by the presence of functional groups and their pair-wise distances to each other make BioPandas a particularly attractive utility library for virtual screening and protein-ligand docking applications.

# References

**McKinney, Wes** (2010). Data Structures for Statistical Computing in Python. *Proceedings of the 9th Python in Science Conference. Vol. 445: pp. 51-56* 

**Berman, Helen, Kim Henrick, and Haruki Nakamura** (2003). Announcing the Worldwide Protein Data Bank. *Nature Structural & Molecular Biology 10.12: pp. 980-980.*

**Berman, H.M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T.N., Weissig, H., Shindyalov, I.N. and Bourne, P.E. (2000)**. The Protein Data Bank. *Nucleic Acids Research 28.1: pp. 235-242.*

**Tripos, L.** (2007). Tripos Mol2 File Format. *St. Louis, MO: Tripos.*
