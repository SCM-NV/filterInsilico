import pubchempy as pcp


def get_data_from_pubchem(identifier, namespace='smiles'):
    """
    """
    df = pcp.get_compounds(identifier, namespace=namespace, as_dataframe=True)

    return df
