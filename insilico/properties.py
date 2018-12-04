from .db import get_data_from_pubchem
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from typing import (Dict, List)
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


def search_property(molecular_properties: List, state: pd.DataFrame,
                    sources: List=["pubchem"]) -> pd.DataFrame:
    """
   Search for a set of `molecular_properties` in different `sources`.
    """
    for key in molecular_properties:
        state[key] = state.smiles.apply(lambda x: search_property_in_sources(x, key, sources))

    return state


def search_property_in_sources(mol: str, prop: str, sources: List):
    """
    Search for molecular `prop` for molecule `mol` in online `sources`.
    """
    funs = {"pubchem": get_data_from_pubchem}
    for s in sources:
        df = funs[s](mol)
        if prop in df:
            return df[prop].values[0]

    return 42


def compute_fingerprint(
        molecules: pd.DataFrame, fingerprint_type: str="topological") -> pd.DataFrame:
    """
    Add a new column to the dataframe with `fingerprint_type`.
    """
    mols = molecules.smiles.apply(lambda x: Chem.MolFromSmiles(x))
    name = 'fingerprint_' + fingerprint_type
    molecules[name] = mols.apply(lambda x: FingerprintMols.FingerprintMol(x))

    return molecules
