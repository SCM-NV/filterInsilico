
from .io import read_molecules
from dask import delayed
from insilico.filters import apply_filter
from insilico.properties import compute_property
from typing import Dict
import pandas as pd


def run_workflow(dict_input: Dict):
    """
    Run the workflow specified in the `dict_input`.
    """
    molecules = read_molecules(dict_input['file_molecules'])
    dag = build_graph(dict_input['steps'], molecules)

    results = runner(dag)

    return results


def build_graph(steps: Dict, molecules: pd.DataFrame) -> object:
    """
    Create a Direct Acyclic Graph containing all the dependencies
    between the filters and the properties to compute.
    """
    # create Dask delayed functions
    delayed_apply_filter = delayed(apply_filter)
    delayed_compute_property = delayed(compute_property)

    for obj in steps:
        if obj.get('apply_filter') is not None:
            print(obj)
        elif obj.get('compute_property') is not None:
            print(obj)

    return 42


def runner(graph: object):
    """
    Run the Direct Acyclic Graph containing all the filters and
    properties.
    """
    pass
