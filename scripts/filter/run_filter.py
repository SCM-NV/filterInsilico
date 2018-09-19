
from filterInsilico.filters import apply_filters
from filterInsilico.process_input import process_input
import argparse


msg = "python run_filter.py -i input.yml"
parser = argparse.ArgumentParser(description=msg)

parser.add_argument('-i', required=True, help="Input file in YAML format")


def main():
    args = parser.parse_args()
    inp = process_input(args.i, 'filter')
    apply_filters(inp)


if __name__ == "__main__":
    main()
