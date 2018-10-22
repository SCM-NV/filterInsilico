from filterInsilico.io import read_molecules
from rdkit import Chem
from typing import Dict
import pandas as pd


def compute_properties(dict_input: Dict, molecules: pd.DataFrame=None) -> pd.DataFrame:
    """
    Calculate a set of molecular properties define in `dict_input`.
    """
    if molecules is None:
        molecules = read_molecules(dict_input['input_file'])

    molecular_properties = dict_input['molecular_properties']

    funs = {'compute_fingerprint': compute_fingerprint}

    for prop in molecular_properties:
        molecules = funs[prop](molecules)

    return molecules


def compute_fingerprint(
        molecules: pd.DataFrame, fingerprint_type: str="topological") -> pd.DataFrame:
    """
    Add a new column to the dataframe with `fingerprint_type`.
    """
    mols = molecules.smiles.apply(lambda x: Chem.MolFromSmiles(x))
    molecules[fingerprint_type] = mols.appy(lambda x: Chem.FingerMols.FingerprintMol(x))

    return molecules
