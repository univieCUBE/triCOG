import pandas as pd
from collections import defaultdict
import sys
import os
import pytest
import shutil
import subprocess
from pathlib import Path

TESTDATA = [
    {"name": "ThreeByThree", "folder": "tests/testdata/ThreeByThree"},
    {"name": "TenPhages", "folder": "tests/testdata/TenPhages"},
    {"name": "Herpes", "folder": "tests/testdata/herpes"}#,
    #{"name": "SevenBySeven", "folder": "tests/testdata/SevenBySeven"}
]

def turn_output_to_dict(output_df):
    dict_output = defaultdict(list)
    cog_ids = output_df["cog_id"].dropna().unique()

    for cog_id in cog_ids:
        filtered = output_df[output_df["cog_id"] == cog_id].sort_values(by=["protein_name", "length", "start", "end"])
        for _, row in filtered.iterrows():
            key = (row["protein_name"], row["length"], row["start"], row["end"])
            dict_output[key].append(filtered.loc[filtered["protein_name"] != row["protein_name"], "protein_name"].values.tolist())

    filtered = output_df[output_df["cog_id"].isna()].sort_values(by=["protein_name", "length", "start", "end"])
    for _, row in filtered.iterrows():
        key = (row["protein_name"], row["length"], row["start"], row["end"])
        dict_output[key].append([])

    for key in dict_output.keys():
        dict_output[key] = [sorted(inner_list) for inner_list in dict_output[key]]
        dict_output[key].sort(key=lambda x: tuple(x) if x else ())

    return dict_output

def run_full_test(folder_path: str, dataset_name: str):
    this_dir = Path(__file__).parent.resolve()
    root_dir = this_dir.parent.resolve() 

    oldCOG_output = pd.read_csv(os.path.join(folder_path, "cls.csv"), header=None, names=["protein_name", "organism", "protein", "length", "start", "end", "cog_id"], index_col=False)

    try:
        testoutput_dir = os.path.join(this_dir, "testdata/tmp_full_run_test")
        shutil.rmtree(testoutput_dir, ignore_errors=True)
        Path(testoutput_dir).mkdir(parents=True, exist_ok=True)

        input_folder = os.path.join(folder_path, "fasta")
        unfiltered_blast = os.path.join(folder_path, "BLASTno")
        filtered_blast = os.path.join(folder_path, "BLASTff")
        evalue_threshold = 0.01
        overlap_threshold = 0.5

        triCOG = os.path.join(root_dir, "triCOG.py")

        result = subprocess.run([
            sys.executable, triCOG,
            "run",
            "-i", input_folder,
            "-u", unfiltered_blast,
            "-f", filtered_blast,
            "-e", str(evalue_threshold),
            "-t", str(overlap_threshold),
            "--orgbyfile",
            "-o", str(testoutput_dir)
        ], capture_output=True, text=True)

        # catch silent CLI failures
        assert result.returncode == 0, f"triCOG.py crashed:\n{result.stderr}"

    except Exception as e:
        raise SystemExit(f"{dataset_name}: An error occurred during the full run v2 test: {e}")

    newCOG_output = pd.read_csv(os.path.join(testoutput_dir, "base_output.csv"), header=0)

    mask = ["protein_name", "length", "start", "end", "cog_id"]

    newCOG_output = newCOG_output[mask]
    oldCOG_output = oldCOG_output[mask]

    newCOG_dict = turn_output_to_dict(newCOG_output)
    oldCOG_dict = turn_output_to_dict(oldCOG_output)

    assert newCOG_dict == oldCOG_dict, f"{dataset_name}: The outputs of newCOG and oldCOG differ."

    shutil.rmtree(testoutput_dir, ignore_errors=True)


@pytest.mark.parametrize("dataset", TESTDATA, ids=[d["name"] for d in TESTDATA])
def test_full_run(dataset):
    """Parametrized test to check the full run with actual processing for different datasets."""
    run_full_test(dataset["folder"], dataset["name"])
    sys.stdout.write(f"{dataset['name']}: passed the full run v2 test.\n")

if __name__ == "__main__":
    for dataset in TESTDATA:
        run_full_test(dataset["folder"], dataset["name"])
        sys.stdout.write(f"{dataset['name']}: passed the full run test.\n")