""" Class for working with stack of PDB files"""

# BioPandas
# Author: Sebastian Raschka <mail@sebastianraschka.com>
# License: BSD 3 clause
# Project Website: http://rasbt.github.io/biopandas/
# Code Repository: https://github.com/rasbt/biopandas
from __future__ import annotations

from biopandas.pdb import PandasPdb
from biopandas.mmcif import PandasMmcif
from biopandas.align import TMAlign
from typing import List, Union, Dict, Callable
import os
from copy import deepcopy


class PandasPdbStack:
    def __init__(self):
        self.pdbs = {}

    def add_pdb(self, source: Union[str, Dict[str, List[str]]], key=None):
        """Adds a single PDB to the stack with automatic or explicit keying.
        :param source: a string which defines a filename, a PDB or UniProt ID or a dictionary with lists of strings.
        :param method: a method to process the source. Auto will automatically determine the source type and processes accordingly.

        :return: None
        """
        if isinstance(source, str):
            # Construct the key from the filename (strip directory and extension)
            if source.endswith(('.pdb', '.ent', '.pdb.gz', '.ent.gz')):
                if os.path.exists(source):
                    self.read_pdb(source, key=key)
                else:
                    raise FileNotFoundError(f"File {source} not found.")
            elif source.endswith(('.cif', 'cif.gz')):
                if os.path.exists(source):
                    self.read_mmcif(source, key=key)
                else:
                    raise FileNotFoundError(f"File {source} not found.")
            elif 'ATOM' in source or 'HETATM' in source:
                self.read_pdb_from_list(source, key)
            elif len(source) in [4, 8] and source.isalnum():  # Simplistic check for PDB ID, also handles new format
                self.fetch_pdb(key, pdb_id=source)
            elif len(source) in [6, 10] and source.isalnum():
                self.fetch_pdb(key, uniprot_id=source)
            else:
                raise ValueError("Unrecognized file format or PDB or UniProt ID.")
        else:
            raise TypeError("Source must be a string or a dictionary with lists of strings.")

    def add_pdbs(self, sources: List[Union[str, List[str], Dict[str, List[str]]]]):
        """Adds multiple PDBs to the BioPandas collection from a list of sources.
        :param sources: a list of multiple PDB files or PDB/UniProt ID-s, that needs to be processed the input types can be mixed.
        :param method: a method to process the sources.

        :return: None
        """
        if not isinstance(sources, dict):
            for source in sources:
                self.add_pdb(source)
        else:
            for key, source in sources.items():
                self.add_pdb(source, key=key)

    def read_pdb(self, file_path: str, key=None):
        """Reads PDB file from disk or URL.
        :param file_path: the path to the PDB file. Reads the file and assigns it to the key of the filename.

        :return: None
        """
        pdb = PandasPdb()
        if key is None:
           key = os.path.splitext(os.path.basename(file_path).replace('.gz', ''))[0]
        pdb.read_pdb(file_path)
        self.pdbs[key] = pdb

    def read_mmcif(self, file_path: str, key=None):
        """Reads mmCIF file from disk or URL.
        :param file_path: the path to the mmCIF file. Reads the file and assigns it to the key of the filename.

        :return: None
        """
        mmcif = PandasMmcif()
        if key is None:
            key = os.path.splitext(os.path.basename(file_path).replace('.gz', ''))[0]
        mmcif.read_mmcif(file_path)
        self.pdbs[key] = mmcif.convert_to_pandas_pdb()

    def fetch_pdb(self, key=None, pdb_id=None, uniprot_id=None):
        """Fetches a PDB file from the RCSB PDB  repository or AF2 database and assigns it to the keyof the same ID.
        :param pdb_id: the PDB ID to fetch.
        :param uniprot_id: the UniProt ID to fetch from the AF2 database.

        :return: None
        """
        pdb = PandasPdb()

        if pdb_id and uniprot_id:
            raise ValueError("Only one of PDB or UniProt ID must be provided.")
        elif pdb_id:
            pdb.fetch_pdb(pdb_id)
            if key is None:
                key = pdb_id
        elif uniprot_id:
            pdb.fetch_pdb(uniprot_id=uniprot_id, source="alphafold2-v4")
            if key is None:
                key = uniprot_id
        else:
            raise ValueError("PDB or UniProt ID must be provided.")
        self.pdbs[key] = pdb

    def read_pdb_from_list(self, pdb_lines: List[str], key):
        """Reads PDB data from a list of lines and assigns it to the specified key.
        :param pdb_lines: a list of PDB lines.
        :param key: a key to associate the PDB with.

        :return: None
        """
        pdb = PandasPdb()
        pdb.read_pdb_from_list(pdb_lines)
        self.pdbs[key] = pdb

    def apply_filter(self, filter_func: Callable, keep_null=True, **kwargs):
        """Applies a filter across all PDBs in the stack and returns a new stack with the filtered PDBs.
        Mandatory inputs for the filter function are the key and the PandasPdb object.
        :param filter_func: a function that processes a pandaspdb objects and returns a modified pandaspdb object.
        :param keep_null: a boolean to modify whether to keep empty structures.

        :return: a new stack with the filtered PDBs.
        """
        new_stack = PandasPdbStack()
        for key, pdb in self.pdbs.items():
            # Apply the filter function which should return a modified PandasPdb object, not just a DataFrame
            new_pdb = deepcopy(pdb)
            filtered_pdb = filter_func(key=key, pdb=new_pdb, **kwargs)  # Ensure that filter_func modifies and returns the entire PandasPdb object
            if keep_null or len(filtered_pdb.df['ATOM']) > 0:
                new_stack.pdbs[key] = filtered_pdb
        return new_stack

    def apply_calculation(self, calculation_func: Callable):
        """Applies a calculation across all PDBs in the stack and returns a dictionary with the calculated values.
        :param calculation_func: a function that processes a pandaspdb objects and returns objects (values, arrays or dictionaries, etc.)

        :return: a dictionary with the calculated values.
        """
        dict_of_values = {}
        for key, pdb in self.pdbs.items():
            # Apply the calculation function which should return a modified PandasPdb object, not just a DataFrame
            value = calculation_func(key, pdb) # Ensure that calculation_func modifies and returns the entire PandasPdb object
            dict_of_values[key] = value
        return dict_of_values

    def delete_entry(self, key):
        """Deletes a PDB entry from the stack.
        :param key: the key of the PDB entry to delete.

        :return: None
        """
        if key in self.pdbs:
            del self.pdbs[key]

    def update_entry(self, key, new_pdb):
        """Updates a PDB entry in the stack.
        :param key: the key of the PDB entry to update.
        :param new_pdb: the new PandasPdb object to associate with the key.

        :return: None
        """
        if key not in self.pdbs:
            self.delete_entry(key)
        self.add_pdb(new_pdb, key=key)


    def write_entries(self, outdir):
        """Writes all PDB entries in the stack to a directory.
        Since there is no pdb to mmcif conversion in biopandas, the output will be in PDB format.
        :param outdir: the directory to write the PDB files to.

        :return: None
        """
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)

        for key, pdb in self.pdbs.items():
            pdb.to_pdb(os.path.join(outdir, f"{key}.pdb"))

    def tmalign_inside(self, target: str=None):
      """For doing TMalign inside a stack, with one of its entries
      :param pdbs: PandasPdbStack with the structures to align. All of them must have only one chain!

      :return: matrix_file_path, tm_score
      """

      # verify that no structure has more than one chain
      for pdb in self.pdbs.values():
          if pdb.df['ATOM']['chain_id'].nunique() > 1:
              raise ValueError("All structures must have only one chain!")

      # if target is provided, check if it is in the stack and use it as the target
      if target:
          if target not in self.pdbs:
              raise ValueError("Target not found in the stack!")
          else:
              target_pdb_id = target
      else:
          # get one structure from the stack - this will be the target. sort by alphabet
          target_pdb_id = sorted(self.pdbs.keys())[0]

      target_pdb = self.pdbs[target_pdb_id]
      target_chain_id = target_pdb.df['ATOM']['chain_id'].unique()[0]

      mobile_pdbs = PandasPdbStack()
      mobile_pdbs.pdbs = {pdb_id: pdb for pdb_id, pdb in self.pdbs.items() if pdb_id != target_pdb_id}

      # get all the chain ID-s in a dictionary
      mobile_chains = {pdb_id: pdb.df['ATOM']['chain_id'].unique()[0] for pdb_id, pdb in self.pdbs.items()}
      tmalign = TMAlign()
      return tmalign.tmalign_to(target_pdb, mobile_pdbs, target_chain_id, mobile_chains)