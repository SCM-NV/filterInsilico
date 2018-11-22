from rdkit import Chem
from typing import (Dict, List)
import pandas as pd


def apply_filter(filters: Dict, molecules: pd.DataFrame=None) -> pd.DataFrame:
    """
    Apply a different set of `filters` to a molecular set.

    :returns: Pandas Dataframe
    """
    keywords = ['functional_groups']
    for key in keywords:
        if key in filters:
            molecules = filter_by_functional_group(molecules, filters[key])

    return molecules


def filter_by_functional_group(molecules: pd.DataFrame, functional_groups: List) -> pd.DataFrame:
    """
    Search for a set of functional_groups
    """
    # Transform molecules to rkdit molecules
    mols = molecules.smiles.apply(lambda x: Chem.MolFromSmiles(x))
    patterns = [Chem.MolFromSmiles(f) for f in functional_groups['smiles']]

    # Check if the functional_groups are in the molecules
    molecules['has_functional_group'] = mols.apply(
        lambda m: any(m.HasSubstructMatch(p) for p in patterns))

    return molecules
