import pandas as pd
import pytest
import os
from COGlse.lse import run_lse_in_memory
from COGprocessinput.read_q2s import read_q2s_to_dict

TEST_DATASETS = [
    {"name": "3By3", "folder": "tests/testdata/ThreeByThree", "description": "Testing ThreeByThree dataset"},
    {"name": "TenPhages", "folder": "tests/testdata/TenPhages", "description": "Starting assertions for TenPhages test case"},
    {"name": "herpes", "folder": "tests/testdata/herpes", "description": "Starting assertions for Herpes test case"}
]

def import_lse_test_data(folder_path: str):
    """Import LSE test data from a folder."""
    expected_lse_file = os.path.join(folder_path, "lse.csv")
    with open(expected_lse_file, "r") as f:
        lse_lines = f.readlines()
        lse_groups = [sorted(line.strip().split(",")) for line in lse_lines]

    prot2org_file = os.path.join(folder_path, "lse_input/orgbyfile_prot2org.csv")
    prot2org_df = pd.read_csv(prot2org_file, header=0)

    hits_file = os.path.join(folder_path, "lse_input/hits.csv")
    hits_df = pd.read_csv(hits_file, header=0)

    q2s_file = os.path.join(folder_path, "lse_input/query2subject.tsv")
    q2s_dict = read_q2s_to_dict(q2s_file) 

    return (lse_groups, prot2org_df, hits_df, q2s_dict)

def run_lse_test(folder_path: str, dataset_name: str):
    """Run LSE test with the provided data."""
    (lse_groups, prot2org_df, hits_df, q2s_dict) = import_lse_test_data(folder_path)
    lse_df = run_lse_in_memory(hits_df, prot2org_df, q2s_dict)

    #sys.stdout.write(f"expected LSE groups: {len(lse_groups)}, actual LSE groups: {len(lse_df)}\n")

    reverse_entity_dict = dict(zip(prot2org_df["protein"], prot2org_df["prot_name"]))

    lse_list = []
    for group in lse_df["group"]:
        group_list = group.split(",")
        group_list = sorted([reverse_entity_dict[int(x)] for x in group_list])
        lse_list.append(group_list)
    
    lse_groups.sort(key=lambda x: x[0])
    lse_list.sort(key=lambda x: x[0])

    assert lse_groups == lse_list, f"LSE groups of {dataset_name} do not match!"

@pytest.mark.parametrize("dataset", TEST_DATASETS, ids=[d["name"] for d in TEST_DATASETS])
def test_lse_match(dataset):
    """Parametrized test to check LSE match for different datasets."""
    run_lse_test(dataset["folder"], dataset["name"])


if __name__ == "__main__":
    for dataset in TEST_DATASETS:
        print(f"\nRunning LSE test for {dataset['name']}: {dataset['description']}")
        run_lse_test(dataset["folder"], dataset["name"])

    