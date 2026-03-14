import sys
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
from COGtriangle.python.data_converter import dataframes_to_dicts


class TestBasicConversion:
    """Test basic DataFrame to dictionary conversions"""
    
    @pytest.fixture
    def basic_dataframes(self):
        """Create basic test dataframes"""
        blast_hits_df = pd.DataFrame({
            'query': [1, 1, 2, 2, 3],
            'subject': [10, 11, 20, 21, 30],
            'evalue': [1e-10, 1e-20, 1e-15, 1e-25, 1e-30],
            'score': [100, 200, 150, 250, 300]
        })
        
        prot2org_df = pd.DataFrame({
            'protein': [1, 2, 3, 10, 11],
            'organism': ['E.coli', 'S.aureus', 'B.subtilis', 'P.aeruginosa', 'K.pneumoniae'],
            'prot_name': ['GeneA', 'GeneB', 'GeneC', 'GeneD', 'GeneE']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1, 2, 3, 10, 11],
            'length': [100, 200, 300, 150, 250]
        })
        
        lse_groups_df = pd.DataFrame({
            'group': ['100,101,102', '200,201', '300']
        })
        
        return blast_hits_df, prot2org_df, prot2len_df, lse_groups_df
    
    def test_blast_hits_conversion(self, basic_dataframes):
        """Test BLAST hits conversion to nested dictionary"""
        blast_df, prot2org_df, prot2len_df, lse_df = basic_dataframes
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        blast_dict = result['blast_hits']
        
        # Check structure
        assert isinstance(blast_dict, dict)
        assert 1 in blast_dict
        assert 2 in blast_dict
        assert 3 in blast_dict
        
        # Check query 1 has 2 hits
        assert len(blast_dict[1]) == 2
        assert blast_dict[1][0]['subject'] == 10
        assert blast_dict[1][1]['subject'] == 11
        
        # Check query 2 has 2 hits
        assert len(blast_dict[2]) == 2
        
        # Check query 3 has 1 hit
        assert len(blast_dict[3]) == 1
        assert blast_dict[3][0]['subject'] == 30
    
    def test_prot2org_conversion(self, basic_dataframes):
        """Test protein to organism mapping"""
        blast_df, prot2org_df, prot2len_df, lse_df = basic_dataframes
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        prot2org = result['prot2org']
        
        assert isinstance(prot2org, dict)
        assert prot2org[1] == 'E.coli'
        assert prot2org[2] == 'S.aureus'
        assert prot2org[3] == 'B.subtilis'
        assert len(prot2org) == 5
    
    def test_prot2len_conversion(self, basic_dataframes):
        """Test protein to length mapping"""
        blast_df, prot2org_df, prot2len_df, lse_df = basic_dataframes
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        prot2len = result['prot2len']
        
        assert isinstance(prot2len, dict)
        assert prot2len[1] == 100
        assert prot2len[2] == 200
        assert prot2len[3] == 300
        assert all(isinstance(k, int) for k in prot2len.keys())
        assert all(isinstance(v, int) for v in prot2len.values())
    
    def test_prot2name_conversion(self, basic_dataframes):
        """Test protein to name mapping"""
        blast_df, prot2org_df, prot2len_df, lse_df = basic_dataframes
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        prot2name = result['prot2name']
        
        assert isinstance(prot2name, dict)
        assert prot2name[1] == 'GeneA'
        assert prot2name[2] == 'GeneB'
        assert prot2name[3] == 'GeneC'
    
    def test_lse_groups_conversion(self, basic_dataframes):
        """Test LSE groups conversion"""
        blast_df, prot2org_df, prot2len_df, lse_df = basic_dataframes
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        all_lse = result['all_lse']
        prot2lse = result['prot2lse']
        
        # Check all_lse structure
        assert isinstance(all_lse, dict)
        assert 100 in all_lse  # First protein of first group
        assert 200 in all_lse  # First protein of second group
        assert 300 in all_lse  # First protein of third group
        
        # Check group contents
        assert all_lse[100] == [100, 101, 102]
        assert all_lse[200] == [200, 201]
        assert all_lse[300] == [300]
        
        # Check prot2lse reverse mapping
        assert prot2lse[100] == 100
        assert prot2lse[101] == 100
        assert prot2lse[102] == 100
        assert prot2lse[200] == 200
        assert prot2lse[201] == 200
        assert prot2lse[300] == 300
    
    def test_return_structure(self, basic_dataframes):
        """Test that all expected keys are in the return dictionary"""
        blast_df, prot2org_df, prot2len_df, lse_df = basic_dataframes
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        expected_keys = {'blast_hits', 'prot2org', 'prot2len', 'prot2name', 'all_lse', 'prot2lse'}
        assert set(result.keys()) == expected_keys


class TestTypeConversions:
    """Test that type conversions work correctly"""
    
    def test_string_protein_ids_converted_to_int(self):
        """Test that string protein IDs are converted to integers"""
        blast_df = pd.DataFrame({
            'query': ['1', '2'],
            'subject': [10, 20]
        })
        
        prot2org_df = pd.DataFrame({
            'protein': ['1', '2'],
            'organism': ['Org1', 'Org2'],
            'prot_name': ['Gene1', 'Gene2']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': ['1', '2'],
            'length': [100, 200]
        })
        
        lse_df = pd.DataFrame({'group': ['1,2']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # Check that keys are integers
        assert all(isinstance(k, int) for k in result['prot2org'].keys())
        assert all(isinstance(k, int) for k in result['prot2len'].keys())
        assert all(isinstance(k, int) for k in result['blast_hits'].keys())
    
    def test_float_lengths_converted_to_int(self):
        """Test that float lengths are converted to integers"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100.0]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert isinstance(result['prot2len'][1], int)
        assert result['prot2len'][1] == 100
    
    def test_organism_names_as_strings(self):
        """Test that organism names are properly converted to strings"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': [123],  # Numeric organism
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert isinstance(result['prot2org'][1], str)
        assert result['prot2org'][1] == '123'


class TestEdgeCases:
    """Test edge cases and boundary conditions"""
    
    def test_empty_blast_hits(self):
        """Test with empty BLAST hits dataframe"""
        blast_df = pd.DataFrame(columns=['query', 'subject', 'evalue'])
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert result['blast_hits'] == {}
    
    def test_single_protein_lse_group(self):
        """Test LSE group with single protein"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['42']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert result['all_lse'][42] == [42]
        assert result['prot2lse'][42] == 42
    
    def test_large_lse_group(self):
        """Test LSE group with many proteins"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        # Create a large group
        large_group = ','.join(str(i) for i in range(1000, 1100))
        lse_df = pd.DataFrame({'group': [large_group]})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert len(result['all_lse'][1000]) == 100
        assert result['prot2lse'][1099] == 1000
    
    def test_multiple_lse_groups(self):
        """Test multiple LSE groups with proper ID assignment"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({
            'group': ['100,101,102', '200,201', '300,301,302,303']
        })
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # Check that first protein is used as LSE ID
        assert 100 in result['all_lse']
        assert 200 in result['all_lse']
        assert 300 in result['all_lse']
        
        # Check that all proteins map back to correct LSE ID
        assert result['prot2lse'][101] == 100
        assert result['prot2lse'][102] == 100
        assert result['prot2lse'][201] == 200
        assert result['prot2lse'][303] == 300
    
    def test_duplicate_proteins_in_different_groups(self):
        """Test behavior when protein appears in multiple LSE groups"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        # Protein 200 appears in both groups
        lse_df = pd.DataFrame({
            'group': ['100,200', '200,201']
        })
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # Last occurrence wins
        assert result['prot2lse'][200] == 200  # From second group
    
    def test_whitespace_in_lse_groups(self):
        """Test LSE groups with whitespace"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        # Groups with spaces (should be handled by int() conversion)
        lse_df = pd.DataFrame({
            'group': [' 100 , 101 , 102 ']
        })
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # int() should strip whitespace
        assert 100 in result['all_lse']
        assert result['all_lse'][100] == [100, 101, 102]
    
    def test_single_hit_per_query(self):
        """Test queries with single BLAST hit"""
        blast_df = pd.DataFrame({
            'query': [1, 2, 3],
            'subject': [10, 20, 30],
            'evalue': [1e-10, 1e-20, 1e-30]
        })
        
        prot2org_df = pd.DataFrame({
            'protein': [1, 2, 3],
            'organism': ['Org1', 'Org2', 'Org3'],
            'prot_name': ['Gene1', 'Gene2', 'Gene3']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1, 2, 3],
            'length': [100, 200, 300]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # Each query should have exactly one hit
        assert all(len(hits) == 1 for hits in result['blast_hits'].values())
    
    def test_many_hits_per_query(self):
        """Test query with many BLAST hits"""
        # Create 100 hits for query 1
        blast_df = pd.DataFrame({
            'query': [1] * 100,
            'subject': list(range(1000, 1100)),
            'evalue': [1e-10] * 100
        })
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert len(result['blast_hits'][1]) == 100


class TestErrorHandling:
    """Test error handling and invalid inputs"""
    
    def test_empty_protein_column_raises_error(self):
        """Test that empty protein column raises ValueError"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [],
            'organism': [],
            'prot_name': []
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        with pytest.raises(ValueError, match="prot2org_df must contain a 'protein' column"):
            dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
    
    def test_missing_columns_raises_error(self):
        """Test that missing required columns raise KeyError"""
        blast_df = pd.DataFrame({'wrong_column': [1]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        with pytest.raises(KeyError):
            dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
    
    def test_invalid_lse_format_raises_error(self):
        """Test that invalid LSE group format raises ValueError"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        # Non-numeric values in group
        lse_df = pd.DataFrame({'group': ['abc,def']})
        
        with pytest.raises(ValueError):
            dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
    
    def test_missing_prot_name_column_raises_error(self):
        """Test that missing prot_name column raises KeyError"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1']
            # Missing prot_name
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        with pytest.raises(KeyError):
            dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)


class TestDataIntegrity:
    """Test data integrity and consistency"""
    
    def test_blast_dict_preserves_all_columns(self):
        """Test that BLAST dict preserves all columns from dataframe"""
        blast_df = pd.DataFrame({
            'query': [1],
            'subject': [10],
            'evalue': [1e-10],
            'bitscore': [100],
            'identity': [95.5],
            'coverage': [100]
        })
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        hit = result['blast_hits'][1][0]
        assert 'subject' in hit
        assert 'evalue' in hit
        assert 'bitscore' in hit
        assert 'identity' in hit
        assert 'coverage' in hit
    
    def test_no_data_loss_in_mappings(self):
        """Test that no data is lost in protein mappings"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1, 2, 3, 4, 5],
            'organism': ['Org1', 'Org2', 'Org3', 'Org4', 'Org5'],
            'prot_name': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1, 2, 3, 4, 5],
            'length': [100, 200, 300, 400, 500]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # All 5 proteins should be in mappings
        assert len(result['prot2org']) == 5
        assert len(result['prot2len']) == 5
        assert len(result['prot2name']) == 5
    
    def test_consistent_protein_ids_across_mappings(self):
        """Test that protein IDs are consistent across all mappings"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        protein_ids = [1, 2, 3, 4, 5]
        prot2org_df = pd.DataFrame({
            'protein': protein_ids,
            'organism': [f'Org{i}' for i in protein_ids],
            'prot_name': [f'Gene{i}' for i in protein_ids]
        })
        
        prot2len_df = pd.DataFrame({
            'protein': protein_ids,
            'length': [i * 100 for i in protein_ids]
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # All mappings should have same protein IDs
        assert set(result['prot2org'].keys()) == set(protein_ids)
        assert set(result['prot2len'].keys()) == set(protein_ids)
        assert set(result['prot2name'].keys()) == set(protein_ids)


class TestRealWorldScenarios:
    """Test realistic scenarios"""
    
    def test_typical_blast_results(self):
        """Test with typical BLAST output structure"""
        # Typical: multiple queries, varying hit counts
        blast_df = pd.DataFrame({
            'query': [1, 1, 1, 2, 2, 3, 4, 4, 4, 4],
            'subject': [10, 11, 12, 20, 21, 30, 40, 41, 42, 43],
            'evalue': [1e-50, 1e-40, 1e-30, 1e-60, 1e-55, 1e-70, 1e-45, 1e-40, 1e-35, 1e-30],
            'pident': [98.5, 95.3, 92.1, 99.2, 97.8, 100.0, 94.5, 93.2, 91.8, 90.5]
        })
        
        prot2org_df = pd.DataFrame({
            'protein': [1, 2, 3, 4, 10, 11, 12, 20, 21, 30, 40, 41, 42, 43],
            'organism': ['E.coli'] * 14,
            'prot_name': [f'Prot{i}' for i in range(14)]
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1, 2, 3, 4, 10, 11, 12, 20, 21, 30, 40, 41, 42, 43],
            'length': [100, 200, 300, 400, 150, 250, 350, 180, 220, 500, 120, 140, 160, 180]
        })
        
        lse_df = pd.DataFrame({
            'group': ['1,2,3,4', '10,11,12', '20,21', '30', '40,41,42,43']
        })
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # Verify structure
        assert len(result['blast_hits']) == 4  # 4 queries
        assert len(result['blast_hits'][1]) == 3  # Query 1 has 3 hits
        assert len(result['blast_hits'][4]) == 4  # Query 4 has 4 hits
        assert len(result['all_lse']) == 5  # 5 LSE groups
    
    def test_cog_clustering_scenario(self):
        """Test scenario similar to COG clustering workflow"""
        # Simulate COG clustering: proteins grouped by orthology
        blast_df = pd.DataFrame({
            'query': [100, 100, 101, 101, 102],
            'subject': [200, 201, 202, 203, 204],
            'evalue': [1e-100] * 5
        })
        
        prot2org_df = pd.DataFrame({
            'protein': [100, 101, 102, 200, 201, 202, 203, 204],
            'organism': ['Org1', 'Org1', 'Org1', 'Org2', 'Org2', 'Org2', 'Org3', 'Org3'],
            'prot_name': [f'COG{i}' for i in range(8)]
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [100, 101, 102, 200, 201, 202, 203, 204],
            'length': [500] * 8  # Similar lengths for orthologs
        })
        
        # LSE groups represent ortholog clusters
        lse_df = pd.DataFrame({
            'group': ['100,200', '101,201,202', '102,203,204']
        })
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        # Verify ortholog grouping
        assert result['prot2lse'][100] == 100
        assert result['prot2lse'][200] == 100  # Same group
        assert result['prot2lse'][101] == 101
        assert result['prot2lse'][202] == 101  # Same group


# Parametrized tests
class TestParametrized:
    """Parametrized tests for comprehensive coverage"""
    
    @pytest.mark.parametrize("protein_ids,expected_count", [
        ([1, 2, 3], 3),
        ([1], 1),
        (list(range(100)), 100),
        ([999, 1000, 1001], 3),
    ])
    def test_various_protein_counts(self, protein_ids, expected_count):
        """Test with various numbers of proteins"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': protein_ids,
            'organism': [f'Org{i}' for i in protein_ids],
            'prot_name': [f'Gene{i}' for i in protein_ids]
        })
        
        prot2len_df = pd.DataFrame({
            'protein': protein_ids,
            'length': [100] * len(protein_ids)
        })
        
        lse_df = pd.DataFrame({'group': ['1']})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert len(result['prot2org']) == expected_count
        assert len(result['prot2len']) == expected_count
    
    @pytest.mark.parametrize("lse_group,expected_lse_id,expected_size", [
        ("100", 100, 1),
        ("100,101", 100, 2),
        ("100,101,102", 100, 3),
        ("1,2,3,4,5,6,7,8,9,10", 1, 10),
    ])
    def test_various_lse_group_sizes(self, lse_group, expected_lse_id, expected_size):
        """Test LSE groups of various sizes"""
        blast_df = pd.DataFrame({'query': [1], 'subject': [10]})
        
        prot2org_df = pd.DataFrame({
            'protein': [1],
            'organism': ['Org1'],
            'prot_name': ['Gene1']
        })
        
        prot2len_df = pd.DataFrame({
            'protein': [1],
            'length': [100]
        })
        
        lse_df = pd.DataFrame({'group': [lse_group]})
        
        result = dataframes_to_dicts(blast_df, prot2org_df, prot2len_df, lse_df)
        
        assert expected_lse_id in result['all_lse']
        assert len(result['all_lse'][expected_lse_id]) == expected_size
