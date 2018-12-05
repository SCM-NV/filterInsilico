from .dependencies import apply_dependencies
from .db import get_data_from_pubchem
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from typing import (Dict, List, Tuple)
import pandas as pd


def compute_property(molecular_properties: Dict, molecules: pd.DataFrame,
                     dependencies: Dict=None) -> pd.DataFrame:
    """
    Calculate a set of `molecular_properties`.

    :param dict molecular_properties: Properties to compute
    :param molecules: Pandas DataFrame containing the properties
    :param dict dependencies: Current task parent
    :returns: Pandas Dataframe

    """
    funs = {'fingerprint': compute_fingerprint}
    if dependencies is not None:
        df = apply_dependencies(molecules, dependencies)
    else:
        df = molecules

    for prop in molecular_properties:
        if prop in funs:
            name, series = funs[prop](df)
            molecules[name] = series

    return molecules


def search_property(molecular_properties: List, molecules: pd.DataFrame,
                    dependencies: Dict=None) -> pd.DataFrame:
    """
    Search for a set of `molecular_properties` in the pubchem database.

    :param dict molecular_properties: Properties to look at
    :param molecules: Pandas DataFrame containing the properties
    :param dict dependencies: Current task parent
    :returns: Pandas Dataframe

    """
    if dependencies is not None:
        df = apply_dependencies(molecules, dependencies)
    else:
        df = molecules

    results = pd.concat(
        {i: search_property_in_sources(s, molecular_properties, "pubchem")
         for i, s in zip(df.index, df.smiles)},
        ignore_index=False)

    results.reset_index(level=1, drop=True, inplace=True)

    for prop in molecular_properties:
        if prop in results.columns:
            molecules[prop] = results[prop]

    return molecules


def search_property_in_sources(mol: str, properties: List, sources: str):
    """
    Search for molecular `properties` for molecule `mol` in online `sources`.
    """
    funs = {"pubchem": get_data_from_pubchem}
    df = funs[sources](mol)

    results = []
    for prop in properties:
        if prop in df:
            results.append(prop)
    return df[results]


def compute_fingerprint(
        molecules: pd.DataFrame, fingerprint_type: str="topological") -> Tuple:
    """
    Add a new column to the dataframe with `fingerprint_type`.
    """
    mols = molecules.smiles.apply(lambda x: Chem.MolFromSmiles(x))
    name = 'fingerprint_' + fingerprint_type
    series = mols.apply(lambda x: FingerprintMols.FingerprintMol(x))
    return name, series
