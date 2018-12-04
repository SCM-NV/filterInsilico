#!/usr/bin/env python

from insilico.runner import run_workflow
from insilico.validate_input import validate_input
import argparse
import pandas as pd

msg = "python run_filter.py -i input.yml"
parser = argparse.ArgumentParser(description=msg)

parser.add_argument('-i', required=True, help="Input file in YAML format")


def main():
    args = parser.parse_args()
    inp = validate_input(args.i, 'filter')
    run_workflow(inp)
    df.as_pickle("results.pkl")


if __name__ == "__main__":
    main()
