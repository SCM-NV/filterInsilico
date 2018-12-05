from typing import Dict
import pandas as pd


def apply_dependencies(state: pd.DataFrame, dependencies: Dict) -> pd.DataFrame:
    """
    """
    if "filters" in dependencies:
        df = state
        for key in dependencies["filters"].keys():
            df = df[df[key]]

        return df
