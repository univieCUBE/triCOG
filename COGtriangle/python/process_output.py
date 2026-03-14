import logging
import pandas as pd

# Import logger
logger = logging.getLogger("cogmaker")

def save_results(results, output_dir, mode, format, filter_size):
    """
    Save results to CSV files in specified mode and output formats
    """
    
    # Build base DataFrame
    cog_entries_data = [
        {
            "protein_name": e.protein_name,
            "organism": e.organism,
            "protein": e.protein,
            "length": e.length,
            "start": e.start,
            "end": e.end,
            "cog_id": e.cog_id
        }
        for e in results.cog_entries
    ]
    cog_entries_data_df = pd.DataFrame(cog_entries_data)
    
    # Switch between modes
    if mode == "single":
        logger.info("Single-COG mode: keeping only largest COG membership per protein")
        base_df = single_COG_mode(cog_entries_data_df)
    elif mode == "multi":
        logger.info("Multi-COG mode: keeping all COG memberships")
        base_df = cog_entries_data_df.copy()
    

    # Filter output by cluster size (number of unique organisms)
    if filter_size != 1:
        logger.info(f"Filtering output to only include clusters with at least {filter_size} unique organisms")
        filtered_df = filter_output(base_df, filter_size)
    else:
        filtered_df = base_df.copy()
   

    # Save output formats
    if "base" in format:
        logger.info("Saving base format")
        filtered_df.to_csv(f"{output_dir}/base_output.csv", index=False)
    if "legacy" in format:
        logger.info("Saving legacy format")
        legacy_df = filtered_df[["protein_name", "organism", "protein_name", "length", "start", "end", "cog_id"]].copy()
        legacy_df.columns = ["protein_name", "organism", "protein_name_copy", "length", "start", "end", "cog_id"]
        # Add comma at the end
        legacy_df[""] = None
        legacy_df.to_csv(f"{output_dir}/legacy_output.csv", index=False, header = False)
    if "minimal" in format:
        logger.info("Saving minimal format")
        minimal_df = filtered_df[["protein_name", "organism", "cog_id"]].copy()
        minimal_df.to_csv(f"{output_dir}/minimal_output.csv", index=False)
    if "multi" in format:
        if mode == "multi":
            multi_df = multi_output(cog_entries_data_df, filter_size)
            multi_df.to_csv(f"{output_dir}/multi_output.csv", index=False)
        else:
            logger.warning("Multi output cannot be generated in single mode.")



def multi_output(base_df: pd.DataFrame, filter_size: int) -> pd.DataFrame:
    """
    Collapses COG results to one line per protein, following the -b logic 
    of the original COGtriangles.reformat.pl script.
    
    A cluster is dropped only if it is a complete subset of exactly one 
    single larger cluster. Proteins in multiple remaining clusters are 
    labeled with a 'multi:' prefix.
    Moreover, there is a filter removing clusters that do not meet the organism
    threshold.
    """
    # Map cluster IDs to set of proteins
    cog_to_prots = base_df.groupby("cog_id")["protein_name"].apply(set).to_dict()
    
    # Rank clusters by num of organisms, num of proteins and cluster ID
    cluster_stats = base_df.groupby("cog_id").agg(
        org_count=('organism', 'nunique'),
        prot_count=('protein_name', 'count')
    )
    
    sorted_ids = cluster_stats.sort_values(
        ['org_count', 'prot_count', 'cog_id'], 
        ascending=[False, False, True]
    ).index.tolist()

    accepted_cogs = []

    # Filter out complete subsets
    for i, cog in enumerate(sorted_ids):
        # Ignore clusters that don't meet the filter size
        if cluster_stats.loc[cog, 'org_count'] < filter_size:
            continue

        current_set = cog_to_prots[cog]
        is_strict_subset = False
        
        # Check against higher ranking clusters
        for previous_cog in accepted_cogs:
            better_set = cog_to_prots[previous_cog]
            if current_set.issubset(better_set):
                is_strict_subset = True
                break
        
        if not is_strict_subset:
            accepted_cogs.append(cog)

    # Reformat output 
    df_filtered = base_df[base_df["cog_id"].isin(accepted_cogs)].copy()
    processed_rows = []
    processed_prot_names = set()

    if not df_filtered.empty:
        # Group by protein to handle multi-assignments
        grouped = df_filtered.groupby("protein_name")
        for prot_name, group in grouped:
            # Sort unique COGs to ensure consistent output strings
            cogs = sorted(group["cog_id"].unique())
            res_row = group.iloc[0].copy()
            
            if len(cogs) > 1:
                # Collapse multiple COGs into a 'multi:' string
                res_row["cog_id"] = "multi:" + ":".join(cogs)
            else:
                res_row["cog_id"] = cogs[0]
                
            processed_rows.append(res_row)
            processed_prot_names.add(prot_name)

    # Handle Singlets
    final_parts = []
    if not pd.DataFrame(processed_rows).empty:
        final_parts.append(pd.DataFrame(processed_rows))

    if filter_size == 1:
        all_prot_names = set(base_df["protein_name"].unique())
        singlet_names = all_prot_names - processed_prot_names
        
        singlet_df = base_df[base_df["protein_name"].isin(singlet_names)].drop_duplicates("protein_name").copy()
        singlet_df["cog_id"] = "" 
        final_parts.append(singlet_df)
    
    if not final_parts:
        return pd.DataFrame(columns=base_df.columns)

    # Merge and Final Sort
    final_df = pd.concat(final_parts, ignore_index=True)
    
    # Final Sort: Clustered proteins first, then singlets, 
    # both ordered by organism and protein name
    clustered_part = final_df[final_df["cog_id"] != ""].sort_values(["organism", "protein_name"])
    singlet_part = final_df[final_df["cog_id"] == ""].sort_values(["organism", "protein_name"])
    
    return pd.concat([clustered_part, singlet_part], ignore_index=True)
    
def single_COG_mode(base_df: pd.DataFrame) -> pd.DataFrame:
    """
    Implement original Perl COGtriangles.reformat.pl strict (-s) mode.
    
    Process
      1. Sort COGs by size (organism count desc, protein count desc, COG ID asc)
      2. Process COGs in that order
      3. Each protein is kept in the first (largest) COG it appears in
      4. COGs that end up too small after deduplication are dropped or demoted
      5. Output: COGs first (sorted by size), then singlets (sorted by organism)
    
    Returns:
        Deduplicated DataFrame with COGs first, singlets last
    """
    if base_df.empty:
        return base_df.copy()
    
    df = base_df.copy()
      
    # Split COGs from singlets
    has_cog = df[df['cog_id'].notna()].copy()
    original_singlets = df[df['cog_id'].isna()].copy()
    
    # If no COGs at all, just return singlets sorted by organism
    if has_cog.empty:
        if not original_singlets.empty:
            return original_singlets.sort_values(['organism', 'protein_name']).reset_index(drop=True)
        return df.copy()
    
    # Calculate original COG statistics
    original_cog_stats = (
        has_cog
        .groupby('cog_id')
        .agg(
            organism_count=('organism', 'nunique'),
            protein_count=('protein_name', 'size')
        )
        .reset_index()
    )
    
    # Sort COGs by size
    original_cog_stats = original_cog_stats.sort_values(
        ['organism_count', 'protein_count', 'cog_id'],
        ascending=[False, False, True]
    )
    
    sorted_cogs = original_cog_stats['cog_id'].tolist()
    
    # Process COGs in order (strict mode logic)
    processed_proteins = set()
    accepted_cog_rows = []
    demoted_singlet_rows = []
    
    for cog_id in sorted_cogs:
        # Get all rows for this COG
        cog_rows = has_cog[has_cog['cog_id'] == cog_id].copy()
        
        # STRICT MODE: Filter out proteins already seen in previous (larger) COGs
        new_rows = cog_rows[~cog_rows['protein_name'].isin(processed_proteins)]
        
        if new_rows.empty:
            # All proteins already assigned to larger COGs - skip this COG
            continue
        
        # Count organisms among the NEW proteins only
        remaining_org_count = new_rows['organism'].nunique()
        
        
        # Check if COG should be demoted to singlets
        # (Perl: if filter_size < 2 and only 1 organism left, demote to singlets)
        if remaining_org_count < 2:
            # Demote: clear cog_id for these proteins
            for idx, row in new_rows.iterrows():
                row_copy = row.copy()
                row_copy['cog_id'] = pd.NA
                demoted_singlet_rows.append(row_copy)
                processed_proteins.add(row['protein_name'])
            continue
        
        # Accept this COG with remaining proteins
        for idx, row in new_rows.iterrows():
            accepted_cog_rows.append(row)
            processed_proteins.add(row['protein_name'])
            
    # Build result parts
    result_parts = []
    
    # Accepted COGs (if any)
    if accepted_cog_rows:
        accepted_df = pd.DataFrame(accepted_cog_rows)
        
        # Re-calculate stats for final COGs (after deduplication)
        final_cog_stats = (
            accepted_df
            .groupby('cog_id')
            .agg(
                organism_count=('organism', 'nunique'),
                protein_count=('protein_name', 'size')
            )
            .reset_index()
        )
        
        # Sort final COGs by size
        final_cog_stats = final_cog_stats.sort_values(
            ['organism_count', 'protein_count', 'cog_id'],
            ascending=[False, False, True]
        )
        
        # Create sort rank
        final_cog_stats['_sort_rank'] = range(len(final_cog_stats))
        
        # Merge and sort
        accepted_df = accepted_df.merge(
            final_cog_stats[['cog_id', '_sort_rank']],
            on='cog_id',
            how='left'
        )
        
        accepted_df = accepted_df.sort_values(
            ['_sort_rank', 'organism', 'protein_name']
        ).drop(columns=['_sort_rank'])
        
        result_parts.append(accepted_df)
    
    # Demoted singlets (if filter_size < 2)
    if demoted_singlet_rows:
        demoted_df = pd.DataFrame(demoted_singlet_rows)
        demoted_df = demoted_df.sort_values(['organism', 'protein_name'])
        result_parts.append(demoted_df)
    
    # Original singlets
    if not original_singlets.empty:
        original_singlets = original_singlets.sort_values(['organism', 'protein_name'])
        result_parts.append(original_singlets)
    
    # Concatenate
    if len(result_parts) == 0:
        return df.iloc[0:0].copy()
    elif len(result_parts) == 1:
        return result_parts[0].reset_index(drop=True)
    else:
        return pd.concat(result_parts, ignore_index=True)

def filter_output(base_df: pd.DataFrame, filter_size: int) -> pd.DataFrame:
    # Count unique organisms per COG
    is_singleton = base_df['cog_id'].isna()
    clusters_df = base_df[~is_singleton].copy()
    if clusters_df.empty:
        return clusters_df
    organism_counts = clusters_df.groupby('cog_id')['organism'].transform('nunique')

    # Keep only clusters that meet filter size
    filtered_df = clusters_df[(organism_counts >= filter_size)].copy()
    
    return filtered_df