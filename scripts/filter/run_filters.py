#!/usr/bin/env python

from filterInsilico.filters import apply_filters
from filterInsilico.properties import compute_properties
from filterInsilico.validate_input import validate_input
import argparse


msg = "python run_filter.py -i input.yml"
parser = argparse.ArgumentParser(description=msg)

parser.add_argument('-i', required=True, help="Input file in YAML format")


def main():
    args = parser.parse_args()
    inp = validate_input(args.i, 'filter')

    # compute_properties
    df = compute_properties(inp)

    # apply filters
    apply_filters(df, inp)


if __name__ == "__main__":
    main()
