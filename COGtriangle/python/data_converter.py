from collections import defaultdict

# Convert dataframe to Python dicts for automatic type conversion to C++
# Could be redone with np.arrays
def dataframes_to_dicts(blast_hits_df, prot2org_df, prot2len_df, lse_groups_df):
    # Raise error if header does not contain 'protein'
    if prot2org_df.empty or "protein" not in prot2org_df.columns or prot2org_df["protein"].empty:
        raise ValueError("prot2org_df must contain a 'protein' column with protein IDs.")

    # Convert hits to list of dictionaries (one per row)
    records = blast_hits_df.to_dict('records')
    # Convert queries to native Python
    query_ids = blast_hits_df['query'].astype(int).tolist()
    # Create BLAST dict, grouped by query ID
    blast_dict = defaultdict(list)
    for i in range(len(records)):
        blast_dict[query_ids[i]].append(records[i])
    
    # Create native Python lists
    p_ids = prot2org_df['protein'].astype(int).tolist()
    p_orgs = prot2org_df['organism'].astype(str).tolist()
    p_names = prot2org_df['prot_name'].tolist()
    
    p_len_ids = prot2len_df['protein'].astype(int).tolist()
    p_lens = prot2len_df['length'].astype(int).tolist()

    # Convert to dictionaries
    prot2org = dict(zip(p_ids, p_orgs))
    prot2len = dict(zip(p_len_ids, p_lens))
    prot2name = dict(zip(p_ids, p_names))

    # Iterate though LSE
    all_lse = {}
    prot2lse = {}
    groups = lse_groups_df["group"].values.tolist()
    for group_str in groups:
        # Check if row is valid string
        if not isinstance(group_str, str): continue
        # Create dictionary with LSE ID as key
        proteins = [int(p) for p in group_str.split(",")]
        lse_id = proteins[0]
        all_lse[lse_id] = proteins
        # Create look-up for each protein
        for prot in proteins:
            prot2lse[prot] = lse_id

    return {
        'blast_hits': dict(blast_dict),
        'prot2org': prot2org,
        'prot2len': prot2len,
        'prot2name': prot2name,
        'all_lse': all_lse,
        'prot2lse': prot2lse
    }