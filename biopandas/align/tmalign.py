import tempfile
import subprocess
import numpy as np
from copy import deepcopy
from biopandas.pdb import PandasPdb
from biopandas.stack.stack import PandasPdbStack


class TMAlign:
    """
    Class to align structures using TMalign and transform the mobile structure(s) while extracting TM-scores.
    TODO: extend to handle multiple chains in multiple structures.
    TODO: modify so that it will be subclass of Align
    """
    def __init__(self, target, mobile, tmalign_path='./USalign'):
        if tmalign_path is None:
            self.tmalign_path = tmalign_path

    def parse_rotation_matrix(self, file_path):
        """Parse the rotation matrix of TMalign and translation vector from the TMalign output file.
        :param file_path: the path to the TMalign output file.

        :return: matrix, translation
        """
        matrix = np.zeros((3, 3))
        translation = np.zeros(3)
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('------ The rotation matrix to rotate'):
                    next(file)  # Skip the header line
                    for i in range(3):
                        parts = next(file).split()
                        translation[i] = float(parts[1])
                        matrix[i, :] = list(map(float, parts[2:5]))
        return matrix, translation

    def filter_and_validate_chain(self, pdb, chain_id):
        """Filter the PandasPdb by chain_id and validate the presence of the chain.
        :param pdb: the PandasPdb object to filter.
        :param chain_id: the chain ID to filter by.

        :return: filtered_pdb
        """
        filtered_pdb = PandasPdb()
        filtered_atoms = pdb.df['ATOM'][pdb.df['ATOM']['chain_id'].isin([chain_id])]
        if filtered_atoms.empty:
            raise ValueError(f"No such chain '{chain_id}' found in the structure.")
        filtered_pdb.df['ATOM'] = filtered_atoms
        return filtered_pdb

    def process_structure_for_tmalign(self, target_file, mobile_pdb, tmalign_path):
        """Handle the TMalign execution and transformation for a given mobile structure and return the transformed mobile structure and TM-score
        :param target_file: the target structure's filepath
        :param mobile_pdb: the mobile structure
        :param tmalign_path: path to tmalign executable

        :return: transformed_mobile, tm_score

        """
        with self.write_pdb_to_temp_file(mobile_pdb) as mobile_file:
            matrix_file_path, tm_score = self.run_tmalign(target_file.name, mobile_file.name, tmalign_path)
            matrix, translation = self.parse_rotation_matrix(matrix_file_path)
            transformed_mobile = deepcopy(mobile_pdb)

            """Apply the rotation matrix and translation vector to the structure."""
            coords = mobile_pdb.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].values
            transformed_coords = np.dot(coords, matrix.T) + translation
            transformed_mobile.df['ATOM'][['x_coord', 'y_coord', 'z_coord']] = transformed_coords
        return transformed_mobile, tm_score

    def write_pdb_to_temp_file(self, pdb):
        """Write a PandasPdb object's data to a temporary PDB file and return the file handle.
        :param pdb: the PandasPdb object to write to the file.

        :return: file handle
        """
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
        pdb.to_pdb(path=temp_file.name, records=None, gz=False, append_newline=True)
        return temp_file

    def tmalign_to(self, target, mobiles,target_chain, mobile_chains):
        """Run TMalign and transform the mobile structure(s) while extracting TM-scores, specifying chains to align.
        :param target: the target structure to align to.
        :param mobiles: the structure(s) to align.
        :param mobile_chains: the chain(s) to align.
        """
        tm_scores = {}
        transformed_structures = {}
        filtered_target_pdb = self.filter_and_validate_chain(target, target_chain)
        with self.write_pdb_to_temp_file(filtered_target_pdb) as target_file:
            if isinstance(mobiles, PandasPdb):
                mobile_atoms = self.filter_and_validate_chain(mobiles, mobile_chains)
                transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, mobile_atoms, self.tmalign_path)
                transformed_structures = transformed_mobile
                tm_scores[mobiles] = tm_score
            elif isinstance(mobiles, PandasPdbStack):
                for key, pdb in mobiles.pdbs.items():
                    selected_chain = mobile_chains[key] if isinstance(mobile_chains, dict) and key in mobile_chains else mobile_chains
                    mobile_atoms = self.filter_and_validate_chain(pdb, selected_chain)
                    transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, pdb, self.tmalign_path)
                    transformed_structures[key] = transformed_mobile
                    tm_scores[key] = tm_score
            else:
                print(type(mobiles))
                raise ValueError("Input must be a PandasPdb object or a PandasPdbStack.")
        return transformed_structures, tm_scores

    def run_tmalign(self, target, mobile):
        """Function to execute TMalign with a rotation matrix output for one target-mobile pair.
        :param target: the structure to align to
        :param mobile: the structure to align
        :param tmalign_path: path to tmalign executable

        :return: matrix_file_path, tm_score
        """
        matrix_file_path = tempfile.mktemp(suffix='.txt')
        command = [self.tmalign_path, target, mobile, '-m', matrix_file_path]
        result = subprocess.run(command, capture_output=True, text=True)

        # Parse the TM-score from stdout
        tm_score = None
        for line in result.stdout.splitlines():
            if line.strip().startswith('TM-score=') and 'Structure_1' in line:
                parts = line.split()
                tm_score = float(parts[1])
                break

        return matrix_file_path, tm_score

    def tmalign_inside(self, pdbs):
      """For doing TMalign inside a stack, with one of its entries
      :param pdbs: PandasPdbStack with the structures to align. All of them must have only one chain!
      :param tmalign_path: path to tmalign executable

      :return: matrix_file_path, tm_score
      """

      # verify that no structure has more than one chain
      for pdb in self.pdbs.values():
          if pdb.df['ATOM']['chain_id'].nunique() > 1:
              raise ValueError("All structures must have only one chain!")

      # get one structure from the stack - this will be the target
      target_pdb_id = list(self.pdbs.keys())[0]

      target_pdb = self.pdbs[target_pdb_id]
      target_chain_id = target_pdb.df['ATOM']['chain_id'].unique()[0]

      mobile_pdbs = PandasPdbStack()
      mobile_pdbs.pdbs = {pdb_id: pdb for pdb_id, pdb in self.pdbs.items() if pdb_id != target_pdb_id}

      # get all the chain ID-s in a dictionary
      mobile_chains = {pdb_id: pdb.df['ATOM']['chain_id'].unique()[0] for pdb_id, pdb in self.pdbs.items()}

      return target_pdb.tmalign_to(mobile_pdbs, target_chain_id, mobile_chains, self.tmalign_path)
