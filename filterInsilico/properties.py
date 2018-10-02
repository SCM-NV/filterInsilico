from typing import Dict
import pandas as pd


def compute_properties(dict_input: Dict) -> pd.DataFrame:
    """
    Calculate a set of molecular properties define in `dict_input`.
    """
    input_file = dict_input['input_file']
    molecules = pd.read_table(input_file, header=None, names=["smiles"])

