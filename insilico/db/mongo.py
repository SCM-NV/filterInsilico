from pymongo import MongoClient
from typing import Dict
import pandas as pd


def store_in_db(db_config: Dict, df: pd.DataFrame) -> object:
    """
    Store a pandas `df` in the database specified in `db_config`.
    """
    db = connect_to_db(db_config['host'], db_config['port'], db_config['db_name'])
    collection = db[db_config['collection_name']]
    return collection.insert_one(df.to_dict()).inserted_id


def connect_to_db(host: str=None, port: int=None, db_name='ligands') -> object:
    """
    Connect to a mongo db instance in `host` using `port` and
    search for  database `name`
    """
    client = MongoClient(host, port)
    db = client[db_name]

    return db
