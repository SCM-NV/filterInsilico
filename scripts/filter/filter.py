from rdkit import Chem
from typing import List
import argparse
import pandas as pd


msg = "filter.py -i <file/containing/the/smiles>"
parser = argparse.ArgumentParser(description=msg)

parser.add_argument('-i', required=True, help="Input file containing the smiles")
parser.add_argument('-f', required)

def main():
    args = parser.parse_args()
    smiles = read_smiles(args.i)
    print(len(list(smiles)))


def read_smiles(file_path: str) -> List:
    """
    Read smiles from file.
    """
    with open(file_path, 'r') as f:
        xs = f.read().split()

        return map(Chem.MolFromSmiles, xs)


if __name__ == "__main__":
    main()
