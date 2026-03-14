
import pandas as pd
from collections import defaultdict

# Read in files
def read_csv(hits,prot2org,prot2len,lse,query2subject=None):
    blast_hits_df = pd.read_csv(hits)
    prot2org_df = pd.read_csv(prot2org)
    prot2len_df = pd.read_csv(prot2len)

    with open(lse, 'r') as f:
        lines = [line.strip().strip('"') for line in f if line.strip()]
    lse_groups_df = pd.DataFrame({'group': lines})

    if query2subject is not None:
        query2subject_map = load_query2subject(query2subject)
    else:
        query2subject_map = {}
    return blast_hits_df, prot2org_df, prot2len_df, lse_groups_df, query2subject_map

# Read query2subject
def load_query2subject(filepath):
    query2subject = defaultdict(set)
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            # Split query from subject(s)
            parts = line.split(maxsplit=1)
            
            if len(parts) < 1:
                continue
                
            try:
                query_id = int(parts[0])
                
                # If there is no subject, skip
                if len(parts) == 1:
                    continue
                
                # Splits subjects 
                subject_parts = parts[1].split(',')
                
                for s_id in subject_parts:
                    s_id = s_id.strip()
                    if s_id:
                        query2subject[query_id].add(int(s_id))
                        
            except ValueError as e:
                print(f"Error in line {line_num}: '{line}' could not be processed ({e})")
                continue
                
    return dict(query2subject)
