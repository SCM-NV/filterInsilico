from typing import Dict
import scipy.io as sio


def read_matlab_db(path_to_db: str) -> Dict:
    """
    Load the data contains in `path_to_db` as a dict.
    """
    return sio.loadmat(path_to_db)
