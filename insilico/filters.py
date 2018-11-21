from filterInsilico.io import read_molecules
from rdkit import Chem
from typing import (Dict, List)
import pandas as pd


def apply_filters(dict_input: Dict, molecules: pd.DataFrame=None) -> pd.DataFrame:
    """
    Apply a different set of filters specified in `dict_input`
    to a molecular set.

    :returns: Pandas Dataframe
    """
    if molecules is None:
        molecules = read_molecules(dict_input['input_file'])
    filters = dict_input['filters']
    if 'functional_groups' in filters:
        df = filter_by_functional_group(molecules, filters['functional_group'])

    return df


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
