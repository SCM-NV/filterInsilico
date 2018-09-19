from rdkit import Chem
from typing import (Dict, List)
import pandas as pd


def apply_filters(dict_input: Dict) -> pd.DataFrame:
    """
    Apply a different set of filters specified in `dict_input`
    to a molecular set.

    :returns: Pandas Dataframe
    """
    input_file = dict_input['input_file']
    molecules = pd.read_table(input_file, header=None, names=["smiles"])

    filters = dict_input['filters']
    if 'functional_groups' in filters:
        df = filter_by_functional_group(molecules, filters['functional_groups'])

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
