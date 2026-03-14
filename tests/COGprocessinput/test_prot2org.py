import pandas as pd
import os
import pytest
from COGprocessinput.map_prot2org import map_prot2org, custom_prot2org_mapping

TEST_DATASETS = [
    {"name": "base example", "folder": "tests/testdata/prot2org/base_example"},
    {"name": "organism by filename", "folder": "tests/testdata/prot2org/orgbyfile"},
    {"name": "no fastas", "folder": "tests/testdata/prot2org/empty"},
    {"name": "missing folder", "folder": "tests/testdata/prot2org/missing_folder"},
    {"name": "custom prot2org mapping", "folder": "tests/testdata/prot2org/custom"}
]

def test_prot2org_basic():
    data = TEST_DATASETS[0]
    p2o_df = map_prot2org(data["folder"])
    expected_p2o = pd.read_csv(os.path.join(data["folder"], "expected_prot2org_basic.csv"), header=0)
    pd.testing.assert_frame_equal(p2o_df, expected_p2o)

def test_prot2org_orgbyfile():
    data = TEST_DATASETS[1]
    p2o_df = map_prot2org(data["folder"], OrgByFile=True)
    expected_p2o = pd.read_csv(os.path.join(data["folder"], "expected_prot2org_orgbyfile.csv"), header=0)
    pd.testing.assert_frame_equal(p2o_df, expected_p2o)

def test_prot2org_no_fastas():
    data = TEST_DATASETS[2]
    with pytest.raises(SystemExit):
        map_prot2org(data["folder"])

def test_prot2org_missing_folder():
    data = TEST_DATASETS[3]
    with pytest.raises(SystemExit):
        map_prot2org(data["folder"])

def test_custom_prot2org_mapping():
    data = TEST_DATASETS[4]
    custom_mapping_path = os.path.join(data["folder"], "custom_p2o.csv")
    p2o_df = custom_prot2org_mapping(custom_mapping_path)
    expected_p2o = pd.read_csv(os.path.join(data["folder"], "expected_custom_p2o.csv"), header=0)
    pd.testing.assert_frame_equal(p2o_df, expected_p2o)

def test_custom_prot2org_mapping_missing_file():
    with pytest.raises(SystemExit):
        custom_prot2org_mapping("non_existent_file.csv")

def test_custom_prot2org_mapping_invalid_format():
    data = TEST_DATASETS[4]
    invalid_format = os.path.join(data["folder"], "invalid_format_1.csv")
    with pytest.raises(SystemExit):
        custom_prot2org_mapping(invalid_format)

def test_custom_prot2org_mapping_empty_file():
    data = TEST_DATASETS[4]
    empty_file = os.path.join(data["folder"], "empty_file.csv")
    with pytest.raises(SystemExit):
        custom_prot2org_mapping(empty_file)