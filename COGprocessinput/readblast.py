import pandas as pd
import os
from typing import cast
import logging

logger = logging.getLogger("readblast")

def find_blast_files(directory:str, suffix:str=".tab", dir_type:str="unfiltered") -> list[str]:
    """Find all BLAST files in the given directory with the specified suffix.

    Args:
        directory (str): directory where the BLAST files are located.
        suffix (str, optional): file suffix to filter BLAST files. Defaults to ".tab".
        dir_type (str, optional): type of directory (unfiltered or filtered). Defaults to "unfiltered".

    Raises:
        SystemExit: if no files with suffix (given by argument "suffix") are found in the given directory.

    Returns:
        list: list of relative paths to BLAST files with the specified suffix.
    """

    # Check if the directory exists
    if not os.path.isdir(directory):
        logger.critical("The specified %s directory %s does not exist.", dir_type, directory)
        raise SystemExit(1)

    # Find all blast files in the given directory with the specified suffix
    blast_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(suffix):
                blast_files.append(os.path.join(root, file))    # find relative path to blast files with specified suffix

    # Check if blast files were found
    if not blast_files:
        logger.critical("No BLAST files with suffix '%s' found in the given %s directory %s.", suffix, dir_type, directory)
        raise SystemExit(1)
    
    return blast_files

def read_and_map_blast_file(blast_files:list[str], entity_dict:dict[str, int], blast_cols:dict[int, str], dir_type:str="unfiltered") -> pd.DataFrame:
    """Read all blast files in a specified directory and map query and subject names (SeqID) to hash integer.

    Args:
        blast_files (list[str]): List of paths to blast files to read and map.
        entity_dict (dict[str, int]): Dictionary mapping original protein names to integer hash values.
        blast_cols (dict[int, str]): Dictionary mapping column indices to column names.
        dir_type (str, optional): Type of directory (unfiltered or filtered). Defaults to "unfiltered".

    Raises:
        SystemExit: If no valid BLAST entries are found after processing.

    Returns:
        pd.DataFrame: DataFrame containing all BLAST hits with mapped query and subject IDs.
    """

    # Initialize list to hold DataFrames for each blast file and total length counter
    hits_df_list = []
    total_length = 0

    # Loop through all blast files
    for blast_file in blast_files:
        logger.info("Reading BLAST file: %s", blast_file)
        
        # Check if the BLAST file is empty
        if os.path.getsize(blast_file) == 0:
            logger.warning("Skipping empty BLAST file: %s", blast_file)
            continue
        
        # Read the BLAST file into a DataFrame; "c" engine is used for faster parsing, and comment lines starting with "#" are ignored
        file_hits = pd.read_csv(blast_file, sep="\t", comment="#", header=None, engine="c")
        
        # Check if the BLAST file has the expected number of columns (BLAST output format 6 or 7 with optional comment lines). If not, skip the file and log a warning.
        if file_hits.shape[1] != 12:
            logger.warning("Skipping BLAST file %s due to unexpected number of columns. Expected %d columns, found %d columns.", blast_file, len(blast_cols), file_hits.shape[1])
            continue

        # Rename columns and map query and subject names to hash integer
        file_hits = file_hits.rename(columns=blast_cols)
        file_hits["query"] = file_hits["query_name"].map(entity_dict)
        file_hits["subject"] = file_hits["subject_name"].map(entity_dict)
        
        # Drop hits with no valid query or subject IDs based on entity mapping
        pre_length = len(file_hits)                                                 # Log the number of hits before dropping invalid entries
        file_hits = file_hits.dropna(subset=["query", "subject"])
        entity_mismatch_count = pre_length - len(file_hits)                         # Calculate the number of hits dropped due to invalid query or subject IDs
        if entity_mismatch_count > 0:
            logger.warning("Dropped %d hits from BLAST file %s due to no valid query or subject IDs based on entity mapping.", entity_mismatch_count, blast_file)
        
        # Keep only relevant columns and ensure correct data types for downstream processing
        file_hits = file_hits[["query", "subject", "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"]]
        file_hits = file_hits.astype({"query": int, "subject": int, "qStart": int, "qEnd": int, "sStart": int, "sEnd": int, "evalue": float, "score": float})
        
        # Append the processed hits to the list and update total length counter
        total_pre_length = total_length
        hits_df_list.append(file_hits)
        total_length += len(file_hits)

        # Warning if no valid BLAST entries were added from the current file after processing
        if total_length == total_pre_length:
            logger.warning("No valid %s BLAST entries were added from file %s after processing.", dir_type, blast_file)

    # Program fails if no valid BLAST entries were found in any of the files after processing, since no BLAST hit means no SymBet triangles
    if not hits_df_list:
        logger.critical("No valid BLAST entries found in %s BLAST files after processing.", dir_type)
        raise SystemExit(1)

    # Concatenate all hits into a single DataFrame
    blast_hits = pd.concat(hits_df_list, ignore_index=True)

    return blast_hits

def process_blast(
        entity_dict:dict[str, int],
        unfiltered_directory:str,
        filtered_directory:str,
        evalue_threshold:float=10.0,
        query_id:int=0,
        subject_id:int=1,
        qStart_id:int=6,
        qEnd_id:int=7,
        sStart_id:int=8,
        sEnd_id:int=9,
        evalue_id:int=10,
        score_id:int=11,
        suffix:str=".tab"
) -> tuple[pd.DataFrame, dict[int, list[int]], pd.DataFrame]:
    """Read and process BLAST files from unfiltered and filtered directories, applying e-value threshold and finding self hits.

    Args:
        entity_dict (dict[str, int]): Dictionary mapping original protein names to integer hash values.
        unfiltered_directory (str): Directory containing unfiltered BLAST files.
        filtered_directory (str): Directory containing filtered BLAST files.
        evalue_threshold (float, optional): E-value threshold to filter BLAST hits. Defaults to 10.0.
        query_id (int, optional): Column index for query name in BLAST output. Defaults to 0.
        subject_id (int, optional): Column index for subject name in BLAST output. Defaults to 1.
        qStart_id (int, optional): Column index for query start position in BLAST output. Defaults to 6.
        qEnd_id (int, optional): Column index for query end position in BLAST output. Defaults to 7.
        sStart_id (int, optional): Column index for subject start position in BLAST output. Defaults to 8.
        sEnd_id (int, optional): Column index for subject end position in BLAST output. Defaults to 9.
        evalue_id (int, optional): Column index for e-value in BLAST output. Defaults to 10.
        score_id (int, optional): Column index for score in BLAST output. Defaults to 11.
        suffix (str, optional): Suffix of BLAST files to be processed. Defaults to ".tab".

    Raises:
        SystemExit: If no valid BLAST entries are found in either unfiltered or filtered BLAST files after processing, or if no self hits can be found for any query and no fallback hits are available.
        SystemExit: If the specified directories do not exist or if no BLAST files with the specified suffix are found in the given directories.

    Returns:
        A tuple containing:
            - pd.DataFrame: The processed unfiltered BLAST hits DataFrame.
            - dict[int, list[int]]: A dictionary mapping query hash values to lists of subject hash values.
            - pd.DataFrame: The processed filtered BLAST hits DataFrame.
    """

    # Define column mapping for BLAST output format 6 or 7
    blast_cols = {query_id: "query_name", subject_id: "subject_name", qStart_id: "qStart", qEnd_id: "qEnd",
                  sStart_id: "sStart", sEnd_id: "sEnd", evalue_id: "evalue", score_id: "score"}


    logger.info("Processing unfiltered BLAST files...")
    
    unfiltered_blast_files = find_blast_files(unfiltered_directory, suffix, dir_type="unfiltered")                      # Find all unfiltered BLAST files with the specified suffix in the given directory
    unfiltered_hits = read_and_map_blast_file(unfiltered_blast_files, entity_dict, blast_cols, dir_type="unfiltered")   # Read and map unfiltered BLAST files to get hits DataFrame with query and subject IDs mapped to hash integers

    # Filter unfiltered BLAST hits based on e-value threshold
    valid_hits = unfiltered_hits[unfiltered_hits["evalue"] <= evalue_threshold].copy()
    skipped_hits = unfiltered_hits[unfiltered_hits["evalue"] > evalue_threshold].copy()

    if len(valid_hits) == 0:    # Check if any valid BLAST entries remain after filtering by e-value threshold; if not, program fails since no BLAST hit means no SymBet triangles
        logger.critical("No valid BLAST entries found in unfiltered BLAST files after applying e-value threshold.")
        raise SystemExit(1)

    logger.debug("Skipped %d hits from unfiltered BLAST files due to e-value threshold.", len(skipped_hits))

    # Process valid unfiltered BLAST hits to create a DataFrame sorted by query and then by score descending; This is one of the three output components and holds the primary SymBets used for triangle construction
    processed_unfiltered_blast = (valid_hits
                                  .sort_values(["query", "score"], ascending=[True, False])
                                  .reset_index(drop=True))



    logger.info("Processing filtered BLAST files...")

    filtered_blast_files = find_blast_files(filtered_directory, suffix, dir_type="filtered")                            # Find all filtered BLAST files with the specified suffix in the given directory
    filtered_hits = read_and_map_blast_file(filtered_blast_files, entity_dict, blast_cols, dir_type="filtered")         # Read and map filtered BLAST files to get hits DataFrame with query and subject IDs mapped to hash integers

    logger.debug("Skipped %d hits from filtered BLAST files due to e-value threshold.", len(filtered_hits[filtered_hits["evalue"] > evalue_threshold]))

    filtered_hits = filtered_hits[filtered_hits["evalue"] <= evalue_threshold]                                          # Filter filtered BLAST hits based on e-value threshold

    # Check if any valid BLAST entries remain in filtered BLAST files after filtering by e-value threshold
    if len(filtered_hits) == 0:
        logger.critical("No valid BLAST entries found in filtered BLAST files after applying e-value threshold.")
        raise SystemExit(1)                                                                                             # SystemExit if no valid BLAST entries are found in filtered BLAST files after processing
    
    # Process filtered BLAST hits to create a dictionary mapping query hash values to lists of subject hash values, sorted by query and then by subject; This is the second of the three output components
    processed_filtered_blast: dict[int, list[int]] = {
        cast(int, query): cast(list[int], subjects)
        for query, subjects in filtered_hits.sort_values(["query", "subject"])
                                            .groupby("query")["subject"]
                                            .apply(list)
                                            .items()
    }



    logger.info("Finding self hits...")
    # Find self hits for each query in the valid unfiltered BLAST hits
    self_hits = valid_hits[valid_hits["query"] == valid_hits["subject"]].copy()
    self_hits["length"] = self_hits["qEnd"] - self_hits["qStart"] + 1                                                   # Calculate length of the hit based on query start and end positions
    self_hits = self_hits[["query", "length", "score"]].rename(columns={"query": "protein"}).reset_index(drop=True)     # Keep only query (renamed to protein), length, and score columns for self hits

    # Find out which queries didn't have self hits
    all_queries = set(unfiltered_hits["query"])
    queries_with_self_hits = set(self_hits["protein"])
    missing_self_hits = all_queries - queries_with_self_hits

    best_fallback_hits = []

    if missing_self_hits:
        valid_hits["length"] = valid_hits["qEnd"] - valid_hits["qStart"] + 1
        # For all queries in valid hits (passed e-value threshold) that have no self hit, find the best hit (highest bit score) as a fallback
        fallback_candidates = (valid_hits[valid_hits["query"].isin(missing_self_hits)]
                               .sort_values(["query", "score"], ascending=[True, False])
                               .drop_duplicates("query"))
        
        if not fallback_candidates.empty:
            # Add the best fallback hits if any are found
            fallback_df = fallback_candidates[["query", "length", "score"]].copy().rename(columns={"query": "protein"})
            best_fallback_hits.append(fallback_df)

            # Update set of queries without a self hit after adding fallback hits
            queries_with_self_hits = queries_with_self_hits | set(fallback_df["protein"])
            missing_self_hits = missing_self_hits - set(fallback_df["protein"])

        # If there are still queries without self hits after adding fallback hits, these queries stem from proteins that have no BLAST hit passing the e-value threshold
        if missing_self_hits and not skipped_hits.empty:
            skipped_hits["length"] = skipped_hits["qEnd"] - skipped_hits["qStart"] + 1
            # Add the best hits from the skipped hits (those that didn't pass the e-value threshold) as fallback for the remaining queries without self hits
            skipped_candidates = (skipped_hits[skipped_hits["query"].isin(missing_self_hits)]
                                    .sort_values(["query", "score"], ascending=[True, False])
                                    .drop_duplicates("query"))
            
            if not skipped_candidates.empty:
                # Add the best fallback hits from the skipped hits if any are found
                skipped_fallback_df = skipped_candidates[["query", "length", "score"]].copy().rename(columns={"query": "protein"})
                best_fallback_hits.append(skipped_fallback_df)

    # Add all fallback hits found to self_hits DataFrame            
    if best_fallback_hits:
        self_hits = pd.concat([self_hits] + best_fallback_hits, ignore_index=True)

    # Sort self hits to ensure consistent order for downstream processing; This is the third output component
    self_hits = (self_hits
                 .sort_values(["protein", "score"], ascending=[True, False])
                 .drop_duplicates(subset=["protein"], keep="first")
                 .reset_index(drop=True))
    
    # Summary for debug mode
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug("Total entries in hits processed: %s", len(processed_unfiltered_blast))
        logger.debug("Unique queries in hits results: %s", len(set(processed_unfiltered_blast["query"])))
        logger.debug("Total entries in query2subject processed: %s", sum([len(v) for v in processed_filtered_blast.values()]))
        logger.debug("Unique queries in query2subject results: %s", len(processed_filtered_blast.keys()))
        logger.debug("Total self hits found: %s", len(self_hits))
        logger.debug("Unique proteins in self hits: %s", len(set(self_hits["protein"])))
    
    return (processed_unfiltered_blast, processed_filtered_blast, self_hits)