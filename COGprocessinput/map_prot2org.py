import os
import pandas as pd
from COGprocessinput.extract_org import extract_organism
import logging

logger = logging.getLogger("prot2org")

def map_prot2org(input_directory:str, OrgByFile:bool = False) -> pd.DataFrame:
    """Create a mapping of protein identifiers to organism names from FASTA files in the specified directory. The organism name is extracted from the FASTA header, or the file name if the -orgbyfile flag is set. Additionally, each protein ID is assigned a unique hashing integer.

    Args:
        input_directory (str): Directory containing the FASTA files on which BLAST was run. All files with endings .fa, .fasta, .fna, .faa will be processed.

    Raises:
        SystemExit: If no FASTA files are found in the specified directory.

    Returns:
        pd.DataFrame: DataFrame containing the protein name (SeqID), identifiers (hash integer) and organism name. The organism name is extracted from the FASTA header, or the file name if the -orgbyfile flag is set. The DataFrame has three columns named: "prot_name", "protein", and "organism".
    """

    logger.info("Mapping proteins to organisms from FASTA files...")

    # Check if the directory exists
    if not os.path.isdir(input_directory):
        logger.critical("The specified input directory %s does not exist.", input_directory)
        raise SystemExit(1)

    input_files = sorted([f for f in os.listdir(input_directory) if f.endswith((".fa", ".fasta", ".fna", ".faa"))])

    if not input_files:
        logger.critical("No FASTA files found in the specified directory: %s", input_directory)
        raise SystemExit(1)
    
    
    # Initialize three lists for each column of the output DataFrame: prot_name (SeqID), protein (hashing integer), and organism
    prot_names, proteins, organisms = [], [], []
    counter = 1
    for file in input_files:
        # Check if the file is empty and skip it if it is
        if os.path.getsize(os.path.join(input_directory, file)) == 0:
            logger.warning("Skipping empty FASTA file: %s", file)
            continue

        file_name = file.rsplit(".", 1)[0]  # Remove file extension to get base name (Only used if --orgbyfile flag is set)
        with open(os.path.join(input_directory, file), 'r') as f:
            for line in f:
                # identify lines starting with ">" and extract the identifier and organism name
                if line[0] == ">":
                    # identifier is the first word after ">" in the header line
                    identifier = line[1:].strip().split()[0]

                    # Organism name is the file name without extension if -orgbyfile flag is set
                    if OrgByFile:
                        organism = file_name
                    else:
                        # Organism name is the last string in square brackets in the header line, if it exists
                        organism = extract_organism(line)
                        
                        if organism == "Unknown":
                            logger.warning("No organism found for identifier %s in file %s. Setting organism to 'Unknown'. This could lead to issues during LSE processing", identifier, file)

                    # Append the identifier, hashing integer, and organism name to the respective lists
                    prot_names.append(identifier)
                    proteins.append(counter)
                    organisms.append(organism)
                    # Increase the counter for the next protein ID
                    counter += 1

    # Create final DataFrame with three columns: prot_name (SeqID), protein (hashing integer), and organism
    prot2org_df = pd.DataFrame({"prot_name": prot_names, "protein": proteins, "organism": organisms})
    return prot2org_df

def custom_prot2org_mapping(custom_prot2org_path:str) -> pd.DataFrame:
    """Create a custom mapping of protein identifiers to organism names from a provided CSV file. The CSV file should have two columns: "protein" (protein ID) and "organism" (organism name).

    Args:
        custom_prot2org_path (str): Path to a custom protein-to-organism mapping file (CSV with columns: protein ID, organism). The CSV is not supposed to have a header.

    Raises:
        SystemExit: If the specified custom mapping file does not exist or cannot be read.

    Returns:
        pd.DataFrame: DataFrame containing the custom protein-to-organism mapping.
    """

    logger.info("Using custom protein-to-organism mapping from file: %s", custom_prot2org_path)

    # Check if the custom mapping file exists
    if not os.path.isfile(custom_prot2org_path):
        logger.critical("The specified custom protein-to-organism mapping file does not exist: %s", custom_prot2org_path)
        raise SystemExit(1)

    try:
        custom_prot2org_df = pd.read_csv(custom_prot2org_path, header=None)
        if custom_prot2org_df.shape[1] != 2:
            logger.critical("The specified custom protein-to-organism mapping file must have exactly two columns: 'protein ID' and 'organism'. Found %d columns.", custom_prot2org_df.shape[1])
            raise SystemExit(1)
        custom_prot2org_df.columns = ["prot_name", "organism"]  # Assign column names
        custom_prot2org_df["protein"] = range(1, len(custom_prot2org_df) + 1)  # Assign unique integer IDs to proteins
        return custom_prot2org_df[["prot_name", "protein", "organism"]]
    
    except pd.errors.EmptyDataError:
        logger.critical("The specified custom protein-to-organism mapping file is empty: %s", custom_prot2org_path)
        raise SystemExit(1)
    except pd.errors.ParserError as e:
        logger.critical("Error parsing custom protein-to-organism mapping file: %s. Please ensure the file is a valid CSV.", e)
        raise SystemExit(1)
    except Exception as e:
        logger.critical("Error reading custom protein-to-organism mapping file into DataFrame using pandas.read_csv(): %s", e)
        raise SystemExit(1)