import pandas as pd


def read_molecules(input_file: str) -> pd.DataFrame:
    return pd.read_table(input_file, header=None, names=["smiles"])
