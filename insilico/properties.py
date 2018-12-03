from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from typing import Dict
import pandas as pd


def compute_property(molecular_properties: Dict, molecules: pd.DataFrame=None) -> pd.DataFrame:
    """
    Calculate a set of molecular properties define in `dict_input`.
    """
    funs = {'fingerprint': compute_fingerprint}

    for prop in molecular_properties:
        if prop in funs:
            molecules = funs[prop](molecules)

    return molecules


def search_property(molecular_properties: Dict, state: pd.DataFrame) -> pd.DataFrame:
    """ """

    return state


def compute_fingerprint(
        molecules: pd.DataFrame, fingerprint_type: str="topological") -> pd.DataFrame:
    """
    Add a new column to the dataframe with `fingerprint_type`.
    """
    mols = molecules.smiles.apply(lambda x: Chem.MolFromSmiles(x))
    name = 'fingerprint_' + fingerprint_type
    molecules[name] = mols.apply(lambda x: FingerprintMols.FingerprintMol(x))

    return molecules
