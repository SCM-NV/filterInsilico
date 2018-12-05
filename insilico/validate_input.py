from .schemas import schema_workflow
from typing import Dict
import yaml


def validate_input(input_file: str) -> Dict:
    """
    Read the `input_file` in YAML format, validate it against the
    corresponding schema and return a nested dictionary with the input.

    :param str input_file: path to the input
    :return: Input as dictionary
    :raise SchemaError: If the input is not valid
    """
    with open(input_file, 'r') as f:
        dict_input = yaml.load(f.read())

    return schema_workflow.validate(dict_input)
