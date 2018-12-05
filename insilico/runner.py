
from .io import read_molecules
from dask import delayed
from insilico.filters import apply_filter
from insilico.properties import (compute_property, search_property)
from typing import (Dict, List)
import dask
import pandas as pd


class DependencyError(Exception):
    pass


def run_workflow(dict_input: Dict) -> pd.DataFrame:
    """
    Run the workflow specified in the `dict_input`.

    :param dict dict_input:  Input provided by the user
    :returns: Pandas DataFrame containing the results
    """
    molecules = read_molecules(dict_input['file_molecules'])
    dag = build_graph(dict_input['steps'], molecules)

    results = runner(dag)

    return results


def build_graph(steps: Dict, state: pd.DataFrame) -> Dict:
    """
    Create a Direct Acyclic Graph containing all the dependencies
    between the filters and the properties to compute.

    :param dict steps: Task to perform
    :param state: Current DataFrame used as state
    :return: Dictionary representing the graph of dependencies
    :raises DependencyError: if the dependencies are incongruent
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
        elif idx in results:
            raise DependencyError("ID for the tasks must be unique")
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


def runner(dag: object) -> pd.DataFrame:
    """
    Run the Direct Acyclic Graph containing all the filters and
    properties.

    :return: Pandas DataFrame containing the results
    """
    results = dask.compute(dag)[0]

    return pd.concat(results, axis=0, sort=False).drop_duplicates().reset_index(drop=True)
