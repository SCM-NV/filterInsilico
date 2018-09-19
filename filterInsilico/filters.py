from rdkit import Chem
from typing import (Dict, List)
import pandas as pd


def apply_filters(dict_input: Dict):
    """
    """
    input_file = dict_input['input_file']
    smiles = read_smiles(input_file)
    print(len(list(smiles)))


def read_smiles(file_path: str) -> List:
    """
    Read smiles from file.
    """
    with open(file_path, 'r') as f:
        xs = f.read().split()

    return map(Chem.MolFromSmiles, xs)
