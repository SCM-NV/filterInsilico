
from .io import read_molecules
from dask import delayed
from insilico.filters import apply_filter
from insilico.properties import (compute_property, search_property)
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

    print(results)
    return results


def build_graph(steps: Dict, state: pd.DataFrame) -> object:
    """
    Create a Direct Acyclic Graph containing all the dependencies
    between the filters and the properties to compute.
    """
    # create Dask delayed functions
    delayed_apply_filter = delayed(apply_filter)
    delayed_compute_property = delayed(compute_property)
    delayed_search_property = delayed(search_property)

    dict_funs = {'apply_filter': delayed_apply_filter,
                 'compute_property': delayed_compute_property,
                 'search_property': delayed_search_property}
    dict_calculations = {
        'apply_filter': 'filters', 'compute_property': 'property', 'search_property': 'property'}

    results = {}
    for obj in steps:
        keywords = ['apply_filter', 'compute_property', 'search_property']
        key = select_calculation(obj, keywords)
        fun = dict_funs[key]
        dict_input = obj[key]
        idx = dict_input['id']
        calc = dict_input[dict_calculations[key]]
        dependencies = dict_input.get('depends_on')
        if not dependencies or dependencies is None:
            results[idx] = fun(calc, state)
        else:
            parent_id = dependencies
            results[idx] = fun(
                calc, results[parent_id], retrieve_dependencies(steps, parent_id, keywords))

    return results


def retrieve_dependencies(steps: dict, parent_id: int, keywords) -> dict:
    """
    Return the dictionary of dependencies of a given task
    """
    def get_dict_task(step):
        for k in keywords:
            if k in step:
                return step[k]

    for step in steps:
        dict_task = get_dict_task(step)
        if dict_task["id"] == parent_id:
            return dict_task


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
    results = dask.compute(dag)[0]

    return pd.concat(results, axis=0, sort=False).drop_duplicates().reset_index(drop=True)
