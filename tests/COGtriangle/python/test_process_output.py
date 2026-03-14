import sys
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
from COGtriangle.python.process_output import single_COG_mode, multi_output


# Fixtures for test data
@pytest.fixture
def empty_df():
    return pd.DataFrame(columns=['protein_name', 'organism', 'cog_id'])

@pytest.fixture
def overlap_data():
    """
    P1 is in two COGs: 
    - C1 (Small: 1 organism)
    - C2 (Large: 2 organisms)
    P2 is only in C2.
    """
    return pd.DataFrame({
        "protein_name": ["P1", "P1", "P2"],
        "organism": ["OrgA", "OrgA", "OrgB"],
        "cog_id": ["C1", "C2", "C2"]
    })

@pytest.fixture
def subset_data():
    """
    COG_B is a complete subset of COG_A.
    """
    return pd.DataFrame({
        "protein_name": ["P1", "P2", "P3", "P1", "P2"],
        "organism": ["Org1", "Org2", "Org3", "Org1", "Org2"],
        "cog_id": ["C_A", "C_A", "C_A", "C_B", "C_B"]
    })


# Tests for mutli_output
class TestMultiOutputMode:

    def test_empty_dataframe(self, empty_df):
        result = multi_output(empty_df, filter_size=1)
        assert result.empty
        assert list(result.columns) == ['protein_name', 'organism', 'cog_id']

    def test_multi_label_for_overlaps(self, overlap_data):
        """In -b mode, P1 should be labeled as multi if both COGs pass filter."""
        # filter_size=1 allows C1 (1 org) to survive
        result = multi_output(overlap_data, filter_size=1)
        
        p1_row = result[result["protein_name"] == "P1"]
        # Subset logic note: C1 is a subset of C2, so C1 is dropped.
        # Thus, P1 just gets C2. 
        # To trigger 'multi:', C1 would need a protein NOT in C2.
        assert p1_row.iloc[0]["cog_id"] == "C2"

    def test_multi_prefix_trigger(self):
        """Test actual multi: prefix when clusters are not strict subsets."""
        df = pd.DataFrame({
            "protein_name": ["P1", "P1", "P2", "P3"],
            "organism": ["OrgA", "OrgA", "OrgB", "OrgC"],
            "cog_id": ["C1", "C2", "C1", "C2"]
        })
        # C1 has {P1, P2}, C2 has {P1, P3}. Neither is a subset of the other.
        result = multi_output(df, filter_size=1)
        p1_row = result[result["protein_name"] == "P1"]
        assert "multi:" in p1_row.iloc[0]["cog_id"]
        assert "C1" in p1_row.iloc[0]["cog_id"]
        assert "C2" in p1_row.iloc[0]["cog_id"]

    def test_filter_size_behavior(self, overlap_data):
        """Clusters below filter_size should not exist or be part of multi strings."""
        # filter_size=2: C1 (1 org) is discarded immediately
        result = multi_output(overlap_data, filter_size=2)
        assert "C1" not in result["cog_id"].to_string()
        assert result[result["protein_name"] == "P1"].iloc[0]["cog_id"] == "C2"


# Tests for single_COG_mode
class TestSingleCOGMode:
    
    def test_empty_dataframe(self, empty_df):
        result = single_COG_mode(empty_df)
        pd.testing.assert_frame_equal(result, empty_df)

    def test_strict_demotion(self, overlap_data):
        """In -s mode, P1 is kept only in the largest COG (C2), C1 is ignored."""
        result = single_COG_mode(overlap_data)
        # P1 and P2 should both be in C2
        assert all(result["cog_id"] == "C2")
        assert len(result) == 2

    def test_subset_logic_comparison(self, subset_data):
        """In both modes, C_B (subset) should effectively disappear."""
        res_single = single_COG_mode(subset_data)
        res_multi = multi_output(subset_data, filter_size=1)
        
        assert "C_B" not in res_single["cog_id"].values
        assert "C_B" not in res_multi["cog_id"].values


# Cross-Mode Comparison
def test_compare_modes_on_overlap():
    """Verify that multi_output combines while single_COG chooses."""
    df = pd.DataFrame({
        "protein_name": ["P1", "P1", "P2", "P3"],
        "organism": ["OrgA", "OrgA", "OrgB", "OrgC"],
        "cog_id": ["C1", "C2", "C1", "C2"]
    })
    
    res_s = single_COG_mode(df)
    res_b = multi_output(df, filter_size=1)
    
    # Strict mode: P1 gets one COG (usually the first/largest found)
    assert ":" not in str(res_s[res_s["protein_name"] == "P1"].iloc[0]["cog_id"])
    
    # Combine mode: P1 gets multi prefix
    assert "multi:" in str(res_b[res_b["protein_name"] == "P1"].iloc[0]["cog_id"])