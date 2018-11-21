from rdkit import Chem
from typing import Dict
import pandas as pd


def compute_property(molecular_properties: Dict, molecules: pd.DataFrame=None) -> pd.DataFrame:
    """
    Calculate a set of molecular properties define in `dict_input`.
    """
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
