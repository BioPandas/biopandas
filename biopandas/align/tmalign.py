import os
import subprocess
import tempfile
from copy import deepcopy

import numpy as np

from biopandas.align.align import Align
from biopandas.pdb import PandasPdb
from biopandas.stack.stack import PandasPdbStack


class TMAlign(Align):
    """
    Class to align structures using TMalign and transform the mobile structure(s) while extracting TM-scores.
    TODO: extend to handle multiple chains in multiple structures.
    """
    def __init__(self, tmalign_path: str=None):
        """
        Initialize the TMAlign object with the path to the TMalign executable.
        :param tmalign_path:

        return None
        """
        super().__init__()

        path_script = os.path.dirname(os.path.abspath(__file__))

        if tmalign_path is not None:
            if os.path.exists(tmalign_path):
                self.tmalign_path = tmalign_path
            else:
                raise FileNotFoundError(f"TMalign executable not found at {tmalign_path}.")
        elif os.path.exists(os.path.join(path_script, './USalign')):
            self.tmalign_path = os.path.join(path_script, './USalign')
        elif os.path.exists(os.path.join(path_script, './USalign.exe')):
            self.tmalign_path = os.path.join(path_script, './USalign.exe')
        # if our script has 'tests' in it, the path changes to "../../biopandas/align/"
        elif os.path.exists(os.path.join(path_script, '../../biopandas/align/USalign')):
            self.tmalign_path = os.path.join(path_script, '../../biopandas/align/USalign')
        elif os.path.exists(os.path.join(path_script, '../../biopandas/align/USalign.exe')):
            self.tmalign_path = os.path.join(path_script, '../../biopandas/align/USalign.exe')
        else:
            raise ValueError("Please provide the path to the TMalign executable.")

    def parse_tmalign_rotation_matrix(self, file_path: str) -> (np.array, np.array):
        """Parse the rotation matrix of TMalign and translation vector from the TMalign output file.
        :param file_path: the path to the TMalign output file.

        :return: matrix, translation
        """
        matrix = np.zeros((3, 3))
        translation = np.zeros(3)
        if os.path.exists(file_path) is False:
            raise FileNotFoundError(f"TMalign output file not found at {file_path}.")

        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('------ The rotation matrix to rotate'):
                    next(file)  # Skip the header line
                    for i in range(3):
                        parts = next(file).split()
                        translation[i] = float(parts[1])
                        matrix[i, :] = list(map(float, parts[2:5]))
        return matrix, translation


    def transform_coords(self, pdb, matrix, translation, type='ATOM'):
        """Apply the rotation matrix and translation vector to the structure.
        :param pdb: the PandasPdb object to transform.
        :param type: the record type to transform.
        :param matrix: the rotation matrix.
        :param translation: the translation vector.

        :return: transformed_pdb
        """
        transformed_pdb = deepcopy(pdb)

        # check if you have x_coord or Cartn_x
        if 'x_coord' in pdb.df[type].columns:
            coord_cols = ['x_coord', 'y_coord', 'z_coord']
        elif 'Cartn_x' in pdb.df[type].columns:
            coord_cols = ['Cartn_x', 'Cartn_y', 'Cartn_z']
        else:
            raise ValueError(f"No recognized coordinate columns found in the {type} dataframe.")

        coords = pdb.df[type][coord_cols].values
        transformed_coords = self.transform(coords, matrix, translation)
        transformed_pdb.df[type][coord_cols] = transformed_coords
        return transformed_pdb


    def process_structure_for_tmalign(self, target_file, mobile_pdb: PandasPdb, mobile_chain: str) -> (PandasPdb, float):
        """Handle the TMalign execution and transformation for a given mobile structure and return the transformed mobile structure and TM-score
        :param target_file: the target structure's filepath
        :param mobile_pdb: the mobile structure

        :return: transformed_mobile, tm_score

        """
        mobile_filtered_pdb = self.filter_and_validate_chain(mobile_pdb, mobile_chain)

        with self.write_pdb_to_temp_file(mobile_filtered_pdb) as mobile_file:
            matrix_file_path, tm_score = self.run_tmalign(target_file.name, mobile_file.name)
            matrix, translation = self.parse_tmalign_rotation_matrix(matrix_file_path)
            transformed_mobile = deepcopy(mobile_pdb)

            """Apply the rotation matrix and translation vector to the structure."""
            transformed_mobile = self.transform_coords(transformed_mobile, type='ATOM', matrix=matrix, translation=translation)
            transformed_mobile = self.transform_coords(transformed_mobile, type='HETATM', matrix=matrix,
                                                       translation=translation)

        # clean up
        os.remove(matrix_file_path) if os.path.exists(matrix_file_path) else None
        os.remove(mobile_file.name) if os.path.exists(mobile_file.name) else None

        return transformed_mobile, tm_score


    def tmalign_to(self, target: PandasPdb,
                   mobiles: [PandasPdb, PandasPdbStack],
                   target_chain: str, mobile_chains: [str, dict]) -> ([PandasPdb, PandasPdbStack], [float, dict]):
        """Run TMalign and transform the mobile structure(s) while extracting TM-scores, specifying chains to align.
        :param target: the target structure to align to, a PandasPdb object.
        :param mobiles: the structure(s) to align, either a PandasPdb object or a PandasPdbStack.
        :param target_chain: the chain of the target structure to align to.
        :param mobile_chains: the chain(s) to align. A dictionary for each structure in the stack or a single chain ID.

        :return: return the transformed structures and the corresponding TM-scores
        """

        filtered_target_pdb = self.filter_and_validate_chain(target, target_chain)
        with self.write_pdb_to_temp_file(filtered_target_pdb) as target_file:
            if isinstance(mobiles, PandasPdb):
                mobile_atoms = self.filter_and_validate_chain(mobiles, mobile_chains)
                transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, mobile_atoms, mobile_chains)
                return transformed_mobile, tm_score
            elif  isinstance(mobiles, PandasPdbStack):
                transformed_stack = PandasPdbStack()
                tm_scores = {}
                for key, mobile_pdb in mobiles.pdbs.items():
                    selected_chain = mobile_chains[key] if isinstance(mobile_chains, dict) and key in mobile_chains else mobile_chains
                    transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, mobile_pdb, selected_chain)
                    transformed_stack.pdbs[key] = transformed_mobile
                    tm_scores[key] = tm_score
                return transformed_stack, tm_scores
            else:
                raise ValueError("Input must be a PandasPdb object or a PandasPdbStack not {type(mobiles)}.")

    def run_tmalign(self, target: str, mobile: str, matrix_file_path: str=None) -> (str, float):
        """Function to execute TMalign with a rotation matrix output for one target-mobile pair.
        :param target: the structure to align to, a filepath
        :param mobile: the structure to align, a filepath

        :return: matrix_file_path, tm_score
        """

        # Verify that the target and mobile structures exist
        if not os.path.exists(target):
            raise FileNotFoundError(f"Target structure not found at {target}.")
        if not os.path.exists(mobile):
            raise FileNotFoundError(f"Mobile structure not found at {mobile}.")

        # If no matrix file path is provided, create a temporary file. Create the directory if it does not exist.
        if matrix_file_path is None:
            matrix_file_path = tempfile.mktemp(suffix='.txt')
        else:
            base_dir = os.path.dirname(matrix_file_path)
            os.makedirs(base_dir, exist_ok=True)

        # Prepare and run the command
        command = [self.tmalign_path, mobile, target, '-m', matrix_file_path]
        result = subprocess.run(command, capture_output=True, text=True)

        # If the process fails, prove
        if result.stderr != '':
            raise ValueError(f"TMalign failed with return code {result.returncode}."
                             f"\nstdout: {result.stdout}"
                             f"\nstderr: {result.stderr}")

        # Parse the TM-score from stdout
        tm_score = None
        for line in result.stdout.splitlines():
            if line.strip().startswith('TM-score=') and 'Structure_1' in line:
                parts = line.split()
                tm_score = float(parts[1])
                break
        return matrix_file_path, tm_score


    def tmalign_in_stack(self, stack: PandasPdbStack, mobile_chains: dict, target: str=None) -> (PandasPdbStack, dict):
        """For doing TMalign inside a stack, with one of its entries
        :param stack: PandasPdbStack with the structures to align. All of them must have only one chain!
        :param target: the target structure to align to. If not provided, the first structure in the stack will be used.

        :return: matrix_file_path, tm_score
        """

        # if target is provided, check if it is in the stack and use it as the target
        if target:
            if target not in stack.pdbs:
                raise ValueError("Target not found in the stack!")
            else:
                target_pdb_id = target
        else:
            # get one structure from the stack - this will be the target. sort by alphabet
            target_pdb_id = sorted(stack.pdbs.keys())[0]

        target_pdb = stack.pdbs[target_pdb_id]
        target_chain_id = target_pdb.df['ATOM']['chain_id'].unique()[0]

        mobile_pdbs = PandasPdbStack()
        mobile_pdbs.pdbs = {pdb_id: pdb for pdb_id, pdb in stack.pdbs.items() if pdb_id != target_pdb_id}

        # align the structures
        transformed_structures, tm_scores = self.tmalign_to(target_pdb, mobile_pdbs, target_chain_id, mobile_chains)

        return target_pdb_id, transformed_structures, tm_scores