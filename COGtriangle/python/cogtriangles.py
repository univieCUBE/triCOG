import sys
from pathlib import Path
import pandas as pd
import logging
from COGtriangle.python.data_converter import dataframes_to_dicts
from COGtriangle.python.read_input import read_csv
from COGtriangle.python.process_output import save_results


# Get path
HERE = Path(__file__).resolve().parent
COG_DIR = HERE.parent  # COGtriangle/
PROJECT_ROOT = COG_DIR.parent  # root/

# Add path
sys.path.insert(0, str(PROJECT_ROOT / "build" / "COGtriangle"))

# Get python logger
logger = logging.getLogger("cogmaker")

"""
Recreates functionality of cogtriangles.cpp
Converts python df to dicts, calls C++ functions
Should be called by main file
"""


def cogtriangles(expect,overlap,blast_hits_df,prot2org_df,prot2len_df,lse_groups_df,query2subject_map,cog_name,start_num,output,mode,format,filter_size,log_level):
    
    # Import pybind11 module and set log level for C++ logger
    import cogmaker
    cogmaker.set_log_level(log_level)
    
    # Convert df to dict
    logger.info("Converting DataFrames...")
    data = dataframes_to_dicts(
        blast_hits_df, prot2org_df, prot2len_df, lse_groups_df
    )
    
    # Debug summary of inputs
    if logger.isEnabledFor(logging.DEBUG):
        summarise_inputs(data,query2subject_map)

    
    # Create HitSet
    logger.info("Creating HitSet...")
    hitset = cogmaker.HitSet(
        expect=expect,
        overlap=overlap,
        blast=data['blast_hits'],
        prot2org=data['prot2org'],
        prot2len=data['prot2len'],
        all_lse=data['all_lse'],
        prot2lse=data['prot2lse'],
        query2subject=query2subject_map,
        cog_name=cog_name,
        start_num=start_num
    )
    
    # Get bidirectional bes hits for LSEs    
    for lse_id, proteins in data['all_lse'].items():
        hitset.insert(set(proteins))

    # Generate COGs
    logger.info("Generating COGs...")
    results = hitset.makeCOGs()

    # Add name back into result
    for entry in results.cog_entries:
        prot_id = int(entry.protein)
        entry.protein_name = data['prot2name'].get(prot_id, "")

    # Save results
    logger.info("Saving results...")
    # Get intermediate edges and add in protein names
    all_edges_data = [
        {
            "from_id": e.from_id,
            "from_name": data['prot2name'].get(int(e.from_id), ""),
            "to_id": e.to_id,
            "to_name": data['prot2name'].get(int(e.to_id), "")
        }
        for e in results.all_edges
    ]
    all_edges_df = pd.DataFrame(all_edges_data)

    # Get intermediate COG edges and add in protein names and COG IDs
    cog_edges_data = [
        {
            "cog_id": e.cog_id,
            "from_name": data['prot2name'].get(int(e.from_id), ""),
            "from_org": data['prot2org'].get(int(e.from_id), ""),
            "to_name": data['prot2name'].get(int(e.to_id), ""),
            "to_org": data['prot2org'].get(int(e.to_id), "")
        }
        for e in results.cog_edges
    ]
    cog_edges_df = pd.DataFrame(cog_edges_data)

    # Save to CSV   
    all_edges_df.to_csv(f"{output}/all_edges.csv", index=False)
    cog_edges_df.to_csv(f"{output}/cog_edges.csv", index=False)
    # Save output
    save_results(results, output, mode, format, filter_size)
    
# Input length for debugging purposes
def summarise_inputs(data, query2subject):
    blastHits = data["blast_hits"]
    prot2org  = data["prot2org"]
    prot2len  = data["prot2len"]
    all_lse   = data["all_lse"]
    prot2lse  = data["prot2lse"]

    logger.debug("SUMMARY: blastHits queries = %d", len(blastHits))
    logger.debug("SUMMARY: prot2org size  = %d", len(prot2org))
    logger.debug("SUMMARY: prot2len size  = %d", len(prot2len))
    logger.debug("SUMMARY: all_lse size   = %d", len(all_lse))
    logger.debug("SUMMARY: prot2lse size  = %d", len(prot2lse))
    logger.debug("SUMMARY: query2subject size = %d", len(query2subject))

    # Log min/max/samples
    if blastHits:
        logger.debug("SUMMARY: blastHits minQuery=%d maxQuery=%d",
                    min(blastHits.keys()), max(blastHits.keys()))
    if all_lse:
        logger.debug("SUMMARY: all_lse minLse=%d maxLse=%d",
                    min(all_lse.keys()), max(all_lse.keys()))


    