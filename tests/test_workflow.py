
from insilico.io import read_molecules
from insilico.properties import (compute_property, search_property)
from insilico.runner import build_graph
import pandas as pd
import yaml

molecules = read_molecules("tests/tests_files/smiles.smi")
path_input = "tests/tests_files/input_test.yml"

with open(path_input, 'r') as f:
    input_file = yaml.load(f.read())


def test_search():
    """
    test interface with pubchem
    """
    df = search_property(['xlogp'], molecules)
    assert isinstance(df, pd.DataFrame)


def test_compute():
    """
    test fingerprint computation
    """
    props = {"property": ["fingerprint"]}
    df = compute_property(props, molecules)

    assert isinstance(df, pd.DataFrame)


def test_graph():
    """
    Check that the graph of dependencies is built correctly.
    """
    dag = build_graph(input_file['steps'], molecules)
    assert isinstance(dag, dict)
