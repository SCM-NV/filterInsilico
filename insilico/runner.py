
from .io import read_molecules
from dask import delayed
from insilico.filters import apply_filter
from insilico.properties import compute_property
from typing import (Dict, List)
import dask
import pandas as pd


def run_workflow(dict_input: Dict):
    """
    Run the workflow specified in the `dict_input`.
    """
    molecules = read_molecules(dict_input['file_molecules'])
    dag = build_graph(dict_input['steps'], molecules)

    results = runner(dag)

    return results


def build_graph(steps: Dict, state: pd.DataFrame) -> object:
    """
    Create a Direct Acyclic Graph containing all the dependencies
    between the filters and the properties to compute.
    """
    # create Dask delayed functions
    delayed_apply_filter = delayed(apply_filter)
    delayed_compute_property = delayed(compute_property)

    dict_funs = {'apply_filter': delayed_apply_filter,
                 'compute_property': delayed_compute_property}
    dict_calculations = {'apply_filter': 'filters',
                         'compute_property': 'property'}

    results = {}
    for obj in steps:
        keyword = select_calculation(obj, ['apply_filter', 'compute_property'])
        fun = dict_funs[keyword]
        dict_input = obj[keyword]
        idx = dict_input['id']
        calc = dict_input[dict_calculations[keyword]]
        dependencies = dict_input.get('depends_on')
        if not dependencies or dependencies is None:
            results[idx] = fun(calc, state)
        else:
            if len(dependencies) == 1:
                parent_id = dependencies[0]
                results[idx] = fun(calc, results[parent_id])
            else:
                raise(NotImplemented)

    return results


def select_calculation(obj: Dict, keywords: List) -> str:
    """
    Select the type of calculation to run.
    """
    for k in keywords:
        if k in obj:
            return k


def runner(dag: object):
    """
    Run the Direct Acyclic Graph containing all the filters and
    properties.
    """
    print(dask.compute(dag)[0])
