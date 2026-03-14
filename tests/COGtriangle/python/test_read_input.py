import sys
from pathlib import Path
import pytest
import pandas as pd
import tempfile
import os
from pathlib import Path
from COGtriangle.python.read_input import read_csv, load_query2subject


class TestLoadQuery2Subject:
    """Tests for the load_query2subject function"""
    
    def test_basic_loading(self, tmp_path):
        """Test basic query2subject file loading"""
        # Create test file
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text(
            "query_id,subject_id\n"
            "1  10,20\n"
            "2  30\n"
        )
        
        result = load_query2subject(str(q2s_file))
        
        assert result == {1: {10, 20}, 2: {30}}
    
    def test_duplicate_mappings(self, tmp_path):
        """Test that duplicate query->subject pairs are handled correctly"""
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text(
            "query_id,subject_id\n"
            "1  10,10\n"
            "1  20\n"
        )
        
        result = load_query2subject(str(q2s_file))
        
        # Sets automatically handle duplicates
        assert result == {1: {10, 20}}
    
    def test_empty_file(self, tmp_path):
        """Test handling of empty file"""
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text("")
        
        result = load_query2subject(str(q2s_file))
        
        assert result == {}
    
    def test_only_header(self, tmp_path):
        """Test file with only header"""
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text("query_id,subject_id\n")
        
        result = load_query2subject(str(q2s_file))
        
        assert result == {}
    
    def test_skip_blank_lines(self, tmp_path):
        """Test that blank lines are skipped"""
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text(
            "query_id,subject_id\n"
            "1  10\n"
            "\n"  # Blank line
            "2  20\n"
            "   \n"  # Whitespace line
        )
        
        result = load_query2subject(str(q2s_file))
        
        assert result == {1: {10}, 2: {20}}
    
    def test_file_not_found(self):
        """Test handling of non-existent file"""
        with pytest.raises(FileNotFoundError):
            load_query2subject("nonexistent_file.csv")


class TestReadCSV:
    """Tests for the read_csv function"""
    
    @pytest.fixture
    def test_files(self, tmp_path):
        """Create a set of valid test files"""
        # BLAST hits
        hits_file = tmp_path / "blast_hits.csv"
        hits_df = pd.DataFrame({
            'query': ['Q1', 'Q2'],
            'subject': ['S1', 'S2'],
            'evalue': [1e-10, 1e-20]
        })
        hits_df.to_csv(hits_file, index=False)
        
        # Protein to organism
        prot2org_file = tmp_path / "prot2org.csv"
        prot2org_df = pd.DataFrame({
            'protein_id': ['P1', 'P2'],
            'organism': ['Org1', 'Org2']
        })
        prot2org_df.to_csv(prot2org_file, index=False)
        
        # Protein to length
        prot2len_file = tmp_path / "prot2len.csv"
        prot2len_df = pd.DataFrame({
            'protein_id': ['P1', 'P2'],
            'length': [100, 200]
        })
        prot2len_df.to_csv(prot2len_file, index=False)
        
        # LSE groups
        lse_file = tmp_path / "lse_groups.txt"
        lse_file.write_text("Group1\nGroup2\nGroup3\n")
        
        # Query to subject mapping
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text(
            "query_id,subject_id\n"
            "1  10\n"
            "2  20\n"
        )
        
        return {
            'hits': str(hits_file),
            'prot2org': str(prot2org_file),
            'prot2len': str(prot2len_file),
            'lse': str(lse_file),
            'query2subject': str(q2s_file)
        }
    
    def test_basic_loading(self, test_files):
        """Test basic file loading without query2subject"""
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            test_files['hits'],
            test_files['prot2org'],
            test_files['prot2len'],
            test_files['lse']
        )
        
        assert len(blast_df) == 2
        assert len(prot2org_df) == 2
        assert len(prot2len_df) == 2
        assert len(lse_df) == 3
        assert q2s_map == {}
    
    def test_with_query2subject(self, test_files):
        """Test loading with query2subject file"""
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            test_files['hits'],
            test_files['prot2org'],
            test_files['prot2len'],
            test_files['lse'],
            test_files['query2subject']
        )
        
        assert q2s_map == {1: {10}, 2: {20}}
    
    def test_lse_strips_whitespace(self, tmp_path, test_files):
        """Test that LSE file strips whitespace and blank lines"""
        lse_file = tmp_path / "lse_with_blanks.txt"
        lse_file.write_text(
            "  Group1  \n"
            "\n"
            "Group2\n"
            "   \n"
            "Group3  \n"
        )
        
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            test_files['hits'],
            test_files['prot2org'],
            test_files['prot2len'],
            str(lse_file)
        )
        
        assert len(lse_df) == 3
        assert list(lse_df['group']) == ['Group1', 'Group2', 'Group3']
    
    def test_empty_lse_file(self, tmp_path, test_files):
        """Test handling of empty LSE file"""
        lse_file = tmp_path / "empty_lse.txt"
        lse_file.write_text("")
        
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            test_files['hits'],
            test_files['prot2org'],
            test_files['prot2len'],
            str(lse_file)
        )
        
        assert len(lse_df) == 0
    
    def test_missing_hits_file(self, test_files):
        """Test error when hits file is missing"""
        with pytest.raises(FileNotFoundError):
            read_csv(
                "nonexistent_hits.csv",
                test_files['prot2org'],
                test_files['prot2len'],
                test_files['lse']
            )
    
    def test_missing_prot2org_file(self, test_files):
        """Test error when prot2org file is missing"""
        with pytest.raises(FileNotFoundError):
            read_csv(
                test_files['hits'],
                "nonexistent_prot2org.csv",
                test_files['prot2len'],
                test_files['lse']
            )
    
    def test_missing_prot2len_file(self, test_files):
        """Test error when prot2len file is missing"""
        with pytest.raises(FileNotFoundError):
            read_csv(
                test_files['hits'],
                test_files['prot2org'],
                "nonexistent_prot2len.csv",
                test_files['lse']
            )
    
    def test_missing_lse_file(self, test_files):
        """Test error when LSE file is missing"""
        with pytest.raises(FileNotFoundError):
            read_csv(
                test_files['hits'],
                test_files['prot2org'],
                test_files['prot2len'],
                "nonexistent_lse.txt"
            )
    
    def test_malformed_csv(self, tmp_path, test_files):
        """Test handling of malformed CSV"""
        bad_csv = tmp_path / "bad_hits.csv"
        bad_csv.write_text("Not,A,Valid\nCSV,,Structure\n")
        
        # pandas will still read this, but the structure may be unexpected
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            str(bad_csv),
            test_files['prot2org'],
            test_files['prot2len'],
            test_files['lse']
        )
        
        # Check that it returns a DataFrame (even if malformed)
        assert isinstance(blast_df, pd.DataFrame)
    
    def test_return_types(self, test_files):
        """Test that all return types are correct"""
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            test_files['hits'],
            test_files['prot2org'],
            test_files['prot2len'],
            test_files['lse'],
            test_files['query2subject']
        )
        
        assert isinstance(blast_df, pd.DataFrame)
        assert isinstance(prot2org_df, pd.DataFrame)
        assert isinstance(prot2len_df, pd.DataFrame)
        assert isinstance(lse_df, pd.DataFrame)
        assert isinstance(q2s_map, dict)


class TestEdgeCases:
    """Test edge cases and boundary conditions"""
    
    def test_query2subject_with_zero_ids(self, tmp_path):
        """Test that zero IDs are handled correctly"""
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text(
            "query_id   subject_id\n"
            "0  0, 1\n"
            "1  0\n"
        )
        
        result = load_query2subject(str(q2s_file))
        
        assert result == {0: {0, 1}, 1: {0}}
    
    def test_query2subject_with_negative_ids(self, tmp_path):
        """Test handling of negative IDs"""
        q2s_file = tmp_path / "query2subject.csv"
        q2s_file.write_text(
            "query_id subject_id\n"
            "-1 10\n"
            "1  -10\n"
        )
        
        result = load_query2subject(str(q2s_file))
        
        # Negative integers are valid
        assert result == {-1: {10}, 1: {-10}}
    
    def test_very_large_ids(self, tmp_path):
        """Test handling of very large integer IDs"""
        q2s_file = tmp_path / "query2subject.csv"
        large_id = 2**31 - 1
        q2s_file.write_text(
            f"query_id  subject_id\n"
            f"{large_id}    {large_id}\n"
        )
        
        result = load_query2subject(str(q2s_file))
        
        assert result == {large_id: {large_id}}
    
    def test_lse_with_only_whitespace_lines(self, tmp_path):
        """Test LSE file with only whitespace"""
        #from read_input import read_csv
        
        lse_file = tmp_path / "whitespace_lse.txt"
        lse_file.write_text("   \n\t\t\n  \n")
        
        # Create minimal other files
        hits_file = tmp_path / "hits.csv"
        pd.DataFrame({'col': [1]}).to_csv(hits_file, index=False)
        
        prot2org_file = tmp_path / "prot2org.csv"
        pd.DataFrame({'col': [1]}).to_csv(prot2org_file, index=False)
        
        prot2len_file = tmp_path / "prot2len.csv"
        pd.DataFrame({'col': [1]}).to_csv(prot2len_file, index=False)
        
        blast_df, prot2org_df, prot2len_df, lse_df, q2s_map = read_csv(
            str(hits_file),
            str(prot2org_file),
            str(prot2len_file),
            str(lse_file)
        )
        
        assert len(lse_df) == 0


# Parametrized tests for comprehensive coverage
class TestParametrized:
    """Parametrized tests for multiple input scenarios"""
    
    @pytest.mark.parametrize("line,expected", [
        ("1 10\n", (1, 10)),
        ("100   200\n", (100, 200)),
        ("0 0\n", (0, 0)),
    ])
    def test_valid_query2subject_lines(self, tmp_path, line, expected):
        """Test various valid line formats"""
        q2s_file = tmp_path / "test.csv"
        q2s_file.write_text(f"query_id,subject_id\n{line}")
        
        result = load_query2subject(str(q2s_file))
        query_id, subject_id = expected
        
        assert query_id in result
        assert subject_id in result[query_id]
    
    
