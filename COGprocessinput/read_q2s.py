import pandas as pd

def read_q2s_to_dict(q2s_path:str) -> dict[int, list[int]]:
    """Reads a query2subject.tsv file (as created by readblast) into the standard data structure used by triCOG.

    Args:
        q2s_path (str): Path to the query2subject.tsv file, which should be a tab-separated file with two columns: the first column contains query IDs (integers) and the second column contains comma-separated lists of subject IDs (integers).

    Returns:
        dict[int, list[int]]: A dictionary mapping each query ID (integer) to a list of subject IDs (integers).
    """
    
    q2s_df = pd.read_csv(q2s_path, sep="\t", header=None)
    q2s_dict = {k: list(map(int, v.split(","))) for k, v in zip(q2s_df.iloc[:, 0], q2s_df.iloc[:, 1])}
    return q2s_dict
