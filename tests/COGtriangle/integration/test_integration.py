import pytest
import subprocess
import shutil
import pandas as pd
from pathlib import Path

# Add used paths
THIS_DIR = Path(__file__).parent.resolve()
DATA_DIR = THIS_DIR / "data"
OUTPUT_DIR = THIS_DIR / "output"
TRITEST_DIR = THIS_DIR.parent
TEST_DIR = TRITEST_DIR.parent
ROOT_DIR = TEST_DIR.parent

@pytest.fixture(scope="session", autouse=True)
def prepare_output_dir():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def get_cluster_map(df):
    """Generates a mapping of cluster ID to member proteins"""
    # Filter unclustered
    clustered = df[df['cluster'].notna() & (df['cluster'] != '')]
    # Group by cluster ID
    groups = clustered.groupby('cluster')['protein'].apply(set).to_list()
    # Use groups for comparison
    return groups

def test_triangles_clustering_matches_original():
    # Remove old output if exists
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Get paths to input
    cmd_hits = "data/hits.csv"
    cmd_prot2org = "data/prot2org.csv"
    cmd_prot2len = "data/prot2len.csv"
    cmd_lse = "data/lse.csv"
    cmd_q2s = "data/query2subject.tsv"
    cmd_output = "output"
    
    
    # Call tool
    cli_script = str(ROOT_DIR / "triCOG.py")

    cmd = [
        "python", cli_script, "triangles",
        "--hits", cmd_hits,
        "--prot2org", cmd_prot2org,
        "--prot2len", cmd_prot2len,
        "--lse", cmd_lse,
        "--q2s", cmd_q2s,
        "-e", "0.01",
        "-t", "0.5",
        "--format", "legacy",
        "--outdir", cmd_output
    ]
    result = subprocess.run(
        cmd, 
        capture_output=True, 
        text=True, 
        cwd=str(THIS_DIR)
    )

    assert result.returncode == 0, f"CLI crashed: {result.stderr}"

    # Paths to outputs
    my_output_csv = THIS_DIR/ "output/legacy_output.csv" 
    original_result = THIS_DIR / "data/original_result.csv"

    # Load results
    names = ['protein', 'organism', 'protname', 'length', 'start', 'end', 'cluster', 'extra']
    df_orig = pd.read_csv(original_result, header=None, names=names)
    df_yours = pd.read_csv(my_output_csv, header=None, names=names)

    
    # Compare results
    orig_clusters = sorted([sorted(list(s)) for s in get_cluster_map(df_orig)])
    your_clusters = sorted([sorted(list(s)) for s in get_cluster_map(df_yours)])

    # Check number of clusters
    assert len(orig_clusters) == len(your_clusters), \
        f"Number of clusters is different! Original: {len(orig_clusters)}, This tool: {len(your_clusters)}"

    # Check cluster content
    if orig_clusters != your_clusters:
        # Print differences
        diff_output = THIS_DIR / "output/comparison_failures.csv"
        # Fin clustering differences
        pytest.fail("Clustering content is different! Check logs for details.")

    # Check unclustered proteins
    orig_unclustered = set(df_orig[df_orig['cluster'].isna() | (df_orig['cluster'] == '')]['protein'])
    your_unclustered = set(df_yours[df_yours['cluster'].isna() | (df_yours['cluster'] == '')]['protein'])
    
    assert orig_unclustered == your_unclustered, \
        f"Differenc in singletons! Missing: {orig_unclustered - your_unclustered}, Extra: {your_unclustered - orig_unclustered}"