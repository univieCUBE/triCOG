import pandas as pd
import sys
import os
import pytest
from COGprocessinput.map_prot2org import map_prot2org
from COGprocessinput.readblast import process_blast

# Test data configurations
TEST_DATASETS = [
    {"name": "3By3", "folder": "tests/testdata/ThreeByThree"},
    {"name": "TenPhages", "folder": "tests/testdata/TenPhages"},
    {"name": "herpes", "folder": "tests/testdata/herpes"},
    {"name": "herpes different suffix", "folder": "tests/testdata/herpes", "suffix": ".blast"}#,
    #{"name": "SevenBySeven", "folder": "tests/testdata/SevenBySeven"}
]

def import_test_data(folder_path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, str, str, str]:
    """Import all expected test data from a folder."""
    # Define file paths
    expected_hits_file = os.path.join(folder_path, "BLASTconv/hits.csv")
    expected_query2subject_file = os.path.join(folder_path, "BLASTconv/query2subject.csv")
    expected_self_file = os.path.join(folder_path, "BLASTconv/self.csv")
    expected_hash_file = os.path.join(folder_path, "BLASTconv/hash.csv")

    # Read expected data into DataFrames
    expected_hits = pd.read_csv(
        expected_hits_file, 
        header=None, 
        names=["query", "subject", "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"],
        dtype={"query": int, "subject": int, "qStart": int, "qEnd": int, 
               "sStart": int, "sEnd": int, "evalue": float, "score": float}
    )
    expected_query2subject = pd.read_csv(
        expected_query2subject_file, 
        header=None, 
        names=["query", "subject"],
        dtype={"query": int, "subject": int}
    )
    expected_self = pd.read_csv(
        expected_self_file, 
        header=None, 
        names=["protein", "length", "score"],
        dtype={"protein": int, "length": int, "score": float}
    )
    expected_hash = pd.read_csv(
        expected_hash_file, 
        header=None, 
        names=["protein", "prot_name"],
        dtype={"protein": int, "prot_name": str}
    )

    # Define paths for BLAST files and input FASTA
    unfiltered_blast_file = os.path.join(folder_path, "BLASTno")
    filtered_blast_file = os.path.join(folder_path, "BLASTff")
    input_fasta_file = os.path.join(folder_path, "fasta")

    return (expected_hits, expected_query2subject, expected_self, expected_hash,
            unfiltered_blast_file, filtered_blast_file, input_fasta_file)


def transform_hits(hits: pd.DataFrame, reverse_entity_dict: dict[int, str], hash_lookup: dict[str, int]) -> pd.DataFrame:
    """Transform hits to match expected format."""
    # Map query and subject IDs back to original protein names and then to hash values of expected output
    hits["query"] = hits["query"].map(reverse_entity_dict).map(hash_lookup)
    hits["subject"] = hits["subject"].map(reverse_entity_dict).map(hash_lookup)
    # sort hits by query, then by evalue ascending, then by score descending so that it can be compared to expected output
    hits = hits.sort_values(by=["query", "evalue", "score"], ascending=[True, True, False]).reset_index(drop=True)
    return hits


def transform_q2s(query2subject: dict[int, list[int]], reverse_entity_dict: dict[int, str], hash_lookup: dict[str, int]) -> pd.DataFrame:
    """Transform query2subject dictionary to DataFrame in expected format."""
    # Read query2subject dict into DataFrame with columns "query" and "subject"
    q2s_list = []
    for query, subjects in query2subject.items():
        for subject in subjects:
            q2s_list.append({"query": query, "subject": subject})
    query2subject_df = pd.DataFrame(q2s_list)

    # Map query and subject IDs back to original protein names and then to hash values of expected output
    query2subject_df["query"] = query2subject_df["query"].map(reverse_entity_dict).map(hash_lookup)
    query2subject_df["subject"] = query2subject_df["subject"].map(reverse_entity_dict).map(hash_lookup)
    # sort query2subject by query and then by subject so that it can be compared to expected output
    query2subject_df = query2subject_df.sort_values(by=["query", "subject"]).reset_index(drop=True)
    return query2subject_df


def transform_self(self: pd.DataFrame, reverse_entity_dict: dict[int, str], hash_lookup: dict[str, int]) -> pd.DataFrame:
    """Transform self DataFrame to match expected format."""
    # Map protein IDs back to original protein names and then to hash values of expected output
    self["protein"] = self["protein"].map(reverse_entity_dict).map(hash_lookup)
    # sort self by protein, then by length, then by score so that it can be compared to expected output
    self = self.sort_values(by=["protein", "length", "score"]).reset_index(drop=True)
    return self


def run_blast_test(dataset_folder: str, dataset_name: str, e_value_threshold: float=0.01, suffix: str=".tab"):
    """
    Generic test runner for process_blast function.

    Args:
        dataset_folder: Path to the test dataset folder
        dataset_name: Name of the dataset for logging purposes
        e_value_threshold: E-value threshold for BLAST processing
    """
    # Import test data
    (expected_hits, expected_q2s, expected_self, expected_hash,
     unfiltered_blast_file, filtered_blast_file, input_fasta_file) = import_test_data(dataset_folder)

    # Process BLAST data
    prot2org_df = map_prot2org(input_fasta_file)
    entity_dict = dict(zip(prot2org_df.iloc[:, 0], prot2org_df.iloc[:, 1]))
    (hits, query2subject, self) = process_blast(
        entity_dict, unfiltered_blast_file, filtered_blast_file, 
        evalue_threshold=e_value_threshold, suffix=suffix
    )

    # Transform results to match expected format
    reverse_entity_dict = dict(zip(prot2org_df.iloc[:, 1], prot2org_df.iloc[:, 0]))
    hash_lookup = dict(zip(expected_hash.iloc[:, 1], expected_hash.iloc[:, 0]))

    hits = transform_hits(hits, reverse_entity_dict, hash_lookup)
    query2subject = transform_q2s(query2subject, reverse_entity_dict, hash_lookup)
    self = transform_self(self, reverse_entity_dict, hash_lookup)

    # Summary
    sys.stdout.write(f"Unique hit queries: {hits['query'].nunique()}\n")
    sys.stdout.write(f"Expected hits queries: {expected_hits['query'].nunique()}\n")
    sys.stdout.write(f"Unique query2subject queries: {query2subject['query'].nunique()}\n")
    sys.stdout.write(f"Expected query2subject queries: {expected_q2s['query'].nunique()}\n")
    sys.stdout.write(f"Unique self proteins: {self['protein'].nunique()}\n")
    sys.stdout.write(f"Expected self proteins: {expected_self['protein'].nunique()}\n")

    # Assertions
    sys.stdout.write(f"{dataset_name}: Starting assertions...\n")

    ## Check lengths
    assert len(hits) == len(expected_hits), f"Expected {len(expected_hits)} hits, got {len(hits)}"
    assert len(query2subject) == len(expected_q2s), f"Expected {len(expected_q2s)} query2subject entries, got {len(query2subject)}"
    assert len(self) == len(expected_self), f"Expected {len(expected_self)} self entries, got {len(self)}"

    sys.stdout.write(
        f"{dataset_name}: Lengths match - {len(hits)} hits, "
        f"{len(query2subject)} query2subject entries, {len(self)} self entries\n"
    )

    ## Check content
    expected_hits = expected_hits.sort_values(
        by=["query", "evalue", "score"], ascending=[True, True, False]
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(hits, expected_hits, check_column_type=False)
    sys.stdout.write(f"{dataset_name}: Hits content matches expected output.\n")

    expected_q2s = expected_q2s.sort_values(by=["query", "subject"]).reset_index(drop=True)
    pd.testing.assert_frame_equal(query2subject, expected_q2s)
    sys.stdout.write(f"{dataset_name}: Query2Subject content matches expected output.\n")

    expected_self = expected_self.sort_values(
        by=["protein", "length", "score"]
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(self, expected_self)
    sys.stdout.write(f"{dataset_name}: Self content matches expected output.\n")


# Parametrized test function for pytest
@pytest.mark.parametrize("dataset", TEST_DATASETS, ids=[d["name"] for d in TEST_DATASETS])
def test_process_blast(dataset: dict[str, str]):
    """Parametrized test for process_blast with different datasets."""
    run_blast_test(dataset["folder"], dataset["name"], suffix=dataset.get("suffix", ".tab"))

def test_empty_directory():
    """Test that process_blast handles an empty directory gracefully."""
    with pytest.raises(SystemExit):
        process_blast({}, "tests/testdata/readblast/empty/BLASTno", "tests/testdata/readblast/empty/BLASTff")

def test_missing_directory():
    """Test that process_blast raises an error when the directory is missing."""
    with pytest.raises(SystemExit):
        process_blast({}, "tests/testdata/readblast/missing_folder/BLASTno", "tests/testdata/readblast/missing_folder/BLASTff")

def test_no_tab_files():
    """Test that process_blast raises an error when no .tab files are present."""
    with pytest.raises(SystemExit):
        process_blast({}, "tests/testdata/readblast/different_suffix/BLASTno", "tests/testdata/readblast/different_suffix/BLASTff")

def test_incorrect_suffix():
    """Test that process_blast raises an error when files have incorrect suffix."""
    with pytest.raises(SystemExit):
        process_blast({}, "tests/testdata/readblast/different_suffix/BLASTno", "tests/testdata/readblast/different_suffix/BLASTff")

if __name__ == "__main__":
    for dataset in TEST_DATASETS:
        sys.stdout.write(f"\nRunning test for dataset: {dataset['name']}\n")
        run_blast_test(dataset["folder"], dataset["name"])