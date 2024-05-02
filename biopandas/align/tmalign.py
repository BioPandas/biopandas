import tempfile
import subprocess
import numpy as np
from copy import deepcopy
from biopandas.pdb import PandasPdb
from biopandas.align.align import Align
import os


class TMAlign(Align):
    """
    Class to align structures using TMalign and transform the mobile structure(s) while extracting TM-scores.
    TODO: extend to handle multiple chains in multiple structures.
    TODO: modify so that it will be subclass of Align
    """
    def __init__(self, tmalign_path=None):
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
        else:
            raise ValueError("Please provide the path to the TMalign executable.")

    def parse_tmalign_rotation_matrix(self, file_path):
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

    def process_structure_for_tmalign(self, target_file, mobile_pdb):
        """Handle the TMalign execution and transformation for a given mobile structure and return the transformed mobile structure and TM-score
        :param target_file: the target structure's filepath
        :param mobile_pdb: the mobile structure

        :return: transformed_mobile, tm_score

        """
        with self.write_pdb_to_temp_file(mobile_pdb) as mobile_file:
            matrix_file_path, tm_score = self.run_tmalign(target_file.name, mobile_file.name)
            matrix, translation = self.parse_tmalign_rotation_matrix(matrix_file_path)
            transformed_mobile = deepcopy(mobile_pdb)

            """Apply the rotation matrix and translation vector to the structure."""
            coords = mobile_pdb.df['ATOM'][['x_coord', 'y_coord', 'z_coord']].values
            transformed_coords = self.transform(coords, matrix, translation)
            transformed_mobile.df['ATOM'][['x_coord', 'y_coord', 'z_coord']] = transformed_coords
        return transformed_mobile, tm_score


    def tmalign_to(self, target:PandasPdb, mobiles: PandasPdb, target_chain: str, mobile_chains: [str, dict]):
        """Run TMalign and transform the mobile structure(s) while extracting TM-scores, specifying chains to align.
        :param target: the target structure to align to, a PandasPdb object.
        :param mobiles: the structure(s) to align, either a PandasPdb object or a PandasPdbStack.
        :param target_chain: the chain of the target structure to align to.
        :param mobile_chains: the chain(s) to align. A dictionary for each structure in the stack or a single chain ID.

        :return: return the transformed structures and the corresponding TM-scores
        """
        tm_scores = {}
        transformed_structures = {}
        filtered_target_pdb = self.filter_and_validate_chain(target, target_chain)
        with self.write_pdb_to_temp_file(filtered_target_pdb) as target_file:
            if isinstance(mobiles, PandasPdb):
                mobile_atoms = self.filter_and_validate_chain(mobiles, mobile_chains)
                transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, mobile_atoms)
                return transformed_mobile, tm_score
            elif type(mobiles).__qualname__ == 'PandasPdbStack':
                for key, pdb in mobiles.pdbs.items():
                    selected_chain = mobile_chains[key] if isinstance(mobile_chains, dict) and key in mobile_chains else mobile_chains
                    mobile_atoms = self.filter_and_validate_chain(pdb, selected_chain)
                    transformed_mobile, tm_score = self.process_structure_for_tmalign(target_file, mobile_atoms)
                    transformed_structures[key] = transformed_mobile
                    tm_scores[key] = tm_score
                return transformed_structures, tm_scores
            else:
                raise ValueError("Input must be a PandasPdb object or a PandasPdbStack not {type(mobiles)}.")

    def run_tmalign(self, target, mobile):
        """Function to execute TMalign with a rotation matrix output for one target-mobile pair.
        :param target: the structure to align to, a filepath
        :param mobile: the structure to align, a filepath

        :return: matrix_file_path, tm_score
        """

        if not os.path.exists(target):
            raise FileNotFoundError(f"Target structure not found at {target}.")
        if not os.path.exists(mobile):
            raise FileNotFoundError(f"Mobile structure not found at {mobile}.")

        matrix_file_path = tempfile.mktemp(suffix='.txt')
        command = [self.tmalign_path, target, mobile, '-m', matrix_file_path]
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


