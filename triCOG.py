import argparse
from pathlib import Path
import faulthandler
import logging
from utils.logging_config import setup_logging

# Get logger
logger = logging.getLogger("COG")
faulthandler.enable()

# Custom type for overlap threshold to ensure it's between 0.0 and 1.0
def overlap_range(value):
    try:
        f_value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")
    
    if f_value < 0.0 or f_value > 1.0:
        raise argparse.ArgumentTypeError(f"{value} is not a valid overlap percentage. Choose in interval [0.0, 1.0].")
    
    return f_value

def main():
    # main parser (if the user calls triCOG.py without subcommand, it will show help)
    parser = argparse.ArgumentParser(
        description=("Triangle-based clustering pipeline\n"
        "\n"
        "Pipeline usage can be found in: python triCOG.py run -h / --help"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Parser for subcommands (processinput, lse, triangles, run)
    subparsers = parser.add_subparsers(dest="command")


   # ------ COGprocessinput ------

    # --- COGprocessinput Subcommand ---
    processinput_parser = subparsers.add_parser(
        "processinput",
        help="Process input FASTA and BLAST files to create intermediate files for LSE and triangle steps"
    )

    input_mutex_group = processinput_parser.add_mutually_exclusive_group(required=True)
    input_mutex_group.add_argument("-i", "--input", help="Input directory containing FASTA files on which the BLAST was run. All FASTA files will be processed.", type=str)
    input_mutex_group.add_argument("--custom-prot2org", help="Custom protein-to-organism mapping file (CSV with columns: protein ID, organism, but no header). If not provided, mapping will be determined from FASTA headers or file names based on --orgbyfile flag.", type=str)

    input_output_group = processinput_parser.add_argument_group("Input/Output", description="Required input and output paths")
    input_output_group.add_argument("-u", "--unfiltered_blast", required=True, help="Directory containing unfiltered BLAST files", type=str)
    input_output_group.add_argument("-f", "--filtered_blast", required=True, help="Directory containing filtered BLAST files", type=str)
    input_output_group.add_argument("-o", "--outdir", help="Output directory for the processed BLAST results", default=".")

    processing_options_group = processinput_parser.add_argument_group("Processing Options", description="Options for processing BLAST files")
    processing_options_group.add_argument("-e", "--evalue_threshold", help="E-value threshold for filtering BLAST results", type=float, default=10.0)
    processing_options_group.add_argument("--orgbyfile", action="store_true", help="Whether to determine organism from file name instead of FASTA header.")
    processing_options_group.add_argument("--suffix", help="File suffix to filter BLAST files", default=".tab")

    processinput_parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging for debugging purposes, output is saved in \"debug.log\"")

    blast_column_group = processinput_parser.add_argument_group("BLAST Column Mapping", description="Specify column indices for relevant fields in BLAST files (0-based, default assumes standard BLAST tabular output format 6)")
    blast_column_group.add_argument("--query_id", help="Column index for query IDs in BLAST files", type=int, default=0, choices=range(12))
    blast_column_group.add_argument("--subject_id", help="Column index for subject IDs in BLAST files", type=int, default=1, choices=range(12))
    blast_column_group.add_argument("--qStart", help="Column index for query start positions in BLAST files", type=int, default=6, choices=range(12))
    blast_column_group.add_argument("--qEnd", help="Column index for query end positions in BLAST files", type=int, default=7, choices=range(12))
    blast_column_group.add_argument("--sStart", help="Column index for subject start positions in BLAST files", type=int, default=8, choices=range(12))
    blast_column_group.add_argument("--sEnd", help="Column index for subject end positions in BLAST files", type=int, default=9, choices=range(12))
    blast_column_group.add_argument("--evalue_id", help="Column index for e-value in BLAST files", type=int, default=10, choices=range(12))
    blast_column_group.add_argument("--score_id", help="Column index for score in BLAST files", type=int, default=11, choices=range(12))

    # --- prot2org.py and readblast.py ---
    def run_processinput(args):
        
        # Clear log file
        with open("debug.log", "w") as f:
            pass
        # Get level
        level = logging.DEBUG if args.verbose else logging.INFO
        setup_logging(level=level)

        # Convert outdir to Path object
        outdir = Path(args.outdir)

        # Warn if directory already exists
        if outdir.exists():
            logger.warning("Warning: Output directory already exists: %s. Files may be overwritten.\n", outdir)

        # Create output directory only after validation
        outdir.mkdir(parents=True, exist_ok=True)

        try:
            # import functions
            from COGprocessinput.map_prot2org import map_prot2org, custom_prot2org_mapping
            from COGprocessinput.readblast import process_blast

            if args.custom_prot2org:
                logger.info("Using custom protein-to-organism mapping from file: %s", args.custom_prot2org)
                prot2org_df = custom_prot2org_mapping(args.custom_prot2org)
                if args.orgbyfile:
                    logger.warning("Warning: --orgbyfile flag is set but will be ignored since a custom protein-to-organism mapping file is provided.")
            else:
                # run prot2org to create mapping dataframe and dictionary
                prot2org_df = map_prot2org(args.input, args.orgbyfile)
            
            entity_dict = dict(zip(prot2org_df.iloc[:, 0], prot2org_df.iloc[:, 1]))

            # run readblast to process BLAST files and create hits dataframe, query2subject mapping, and self hits dataframe
            (processed_unfiltered_blast, processed_filtered_blast, self_hits) = process_blast(
                entity_dict=entity_dict,
                unfiltered_directory=args.unfiltered_blast,
                filtered_directory=args.filtered_blast,
                query_id=args.query_id,
                qStart_id=args.qStart,
                qEnd_id=args.qEnd,
                subject_id=args.subject_id,
                sStart_id=args.sStart,
                sEnd_id=args.sEnd,
                evalue_id=args.evalue_id,
                evalue_threshold=args.evalue_threshold,
                score_id=args.score_id,
                suffix=args.suffix
            )

            # write output to output directory
            processed_unfiltered_blast.to_csv(f"{args.outdir}/hits.csv", index=False)

            with open(f"{args.outdir}/query2subject.tsv", "w") as f:
                for query, targets in processed_filtered_blast.items():
                    f.write(f"{query}\t" + ",".join(map(str, targets)) + "\n")
            
            self_hits.to_csv(f"{args.outdir}/self_hits.csv", index=False)     
            
            prot2org_df.to_csv(f"{args.outdir}/prot2org.csv", index=False, header=["prot_name", "protein", "organism"])

            logger.info("COGprocessinput completed successfully. Output written to %s", args.outdir)
        
        except Exception as e:
            logger.error("Error while running COGprocessinput: %s", e)

    processinput_parser.set_defaults(func=run_processinput)

    # ------ lse ------

    # --- lse Subcommand ---
    lse_parser = subparsers.add_parser(
        "lse",
        help="Build lineage-specific expansion groups"
    )
    lse_parser.add_argument("--hits", required=True, help="Input hits file (hits.csv / columns: query,subject,evalue,score)")
    lse_parser.add_argument("--prot2org", required=True, help="Protein-to-organism mapping CSV (must contain columns: protein,organism).")
    lse_parser.add_argument("--job", default=None, help="Optional job file (query_org,outgroup_org). If omitted, default all-vs-all outgroups are used.",)
    lse_parser.add_argument("--q2s", required=True, help="Query-to-subject .tsv (query<TAB>subject1,subject2,...). Used for strict presence filtering.")   
    lse_parser.add_argument("--output", default="lse.csv", help="Output file (one group per line: id1,id2,...) [default: lse.csv]")
    lse_parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging for debugging purposes, output is saved in \"debug.log\"")


    def run_lse_cmd(args):
        # Clear log file
        with open("debug.log", "w") as f:
            pass
        # Get level
        level = logging.DEBUG if args.verbose else logging.INFO
        setup_logging(level=level)

        try:
            from COGlse import lse

            df = lse.run_lse(
                hits_path=args.hits,
                prot2org_path=args.prot2org,
                job_path=args.job,
                presence_map=args.q2s,
                output_path=args.output
                )

            print(f"LSE completed successfully. Groups written to: {args.output}  (n_groups={len(df)})")
        except Exception as e:
            print("Error in lse: %s" % e)

    lse_parser.set_defaults(func=run_lse_cmd)



    # ------ triangles ------

    # --- triangles Subcommand ---
    triangles_parser = subparsers.add_parser(
        "triangles",
        help="Build clusters using triangle algorithm",
        formatter_class=argparse.RawTextHelpFormatter
    )
    triangles_parser.add_argument("--hits", required=True, help="Input hits file (processed_unfiltered_blast.csv)")
    triangles_parser.add_argument("--prot2org", required=True, help="Protein-to-organism mapping file (name,ID,organism)")
    triangles_parser.add_argument("--prot2len",required=True, help="Protein-to-length mapping file (name, length)")
    triangles_parser.add_argument("--lse", required=True,help="Lineage-specific expansion file (one line per LSE)")
    triangles_parser.add_argument("--q2s", required=False,help="Filter file to limit possible edges (query, subjects)")
    triangles_parser.add_argument("-e", "--evalue", type=float, default=10.0, help="E-value threshold for BLAST filtering (default: 10.0)")
    triangles_parser.add_argument("-t", "--overlap", type=overlap_range, default=0.75, help="Overlap threshold (default: 0.75)")
    triangles_parser.add_argument("-n", "--cog_name", type=str, default="COG", help="Prefix for cluster names (default: COG)")
    triangles_parser.add_argument("-s", "--start_number", type=int, default=1, help="Starting number of cluster names (default: 1)")
    triangles_parser.add_argument("-o", "--outdir",required=True, help="Output directory")
    triangles_parser.add_argument("--mode", choices=["multi","single"], default="multi", help="Allow a protein to be in more than one COG (multi) or keep only occurrence in biggest COG (single)")
    triangles_parser.add_argument("--format", choices=["legacy","base","minimal","multi"],default="base",nargs="+", help="Output formats:\n" "base = protein name, organism, protein ID, length, start, end, cluster ID;\n""legacy = protein name, organism, protein name, length, start, end, cluster ID;\n" "minimal = protein name, organism, cluster ID;\n" "multi = Collapsed form with only one line per protein and \"multi\" keyword for proteins in more than one cluster")
    triangles_parser.add_argument("-fs","--filter_size", type=int, default=1, help="Minimum number of unique organisms required for a cluster to be included in the output.\nIf not set, no filtering is applied.")
    triangles_parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging for debugging purposes, output is saved in \"debug.log\"")

    def run_triangles_cmd(args):
        from COGtriangle.python.cogtriangles import cogtriangles
        from COGtriangle.python.read_input import read_csv

        # Clear log file
        with open("debug.log", "w") as f:
            pass
        # Get level
        level = logging.DEBUG if args.verbose else logging.INFO
        setup_logging(level=level)

        try: 
            # Load DataFrames from arg file paths (q2s optional)
            blast_hits_df, prot2org_df, prot2len_df, lse_groups_df, query2subject_map = read_csv(
                args.hits, args.prot2org, args.prot2len, args.lse, args.q2s
            )
            # Run triangle algorithm
            df = cogtriangles(args.evalue, args.overlap, blast_hits_df, prot2org_df, prot2len_df, lse_groups_df, query2subject_map, args.cog_name, args.start_number, args.outdir, args.mode, args.format, args.filter_size, level)
            logger.info("Triangle algorithm run successfully")

        except Exception as e:
            logger.error("Error in triangles: %s", e)
    triangles_parser.set_defaults(func=run_triangles_cmd)




    # ------ run ------

    # --- run Subcommand (Pipeline: processinput → lse → triangles) ---
    run_parser = subparsers.add_parser(
        "run",
        help="Run pipeline: processinput → lse → triangles",
        formatter_class=argparse.RawTextHelpFormatter
    )

    input_mutex_group = run_parser.add_mutually_exclusive_group(required=True)
    input_mutex_group.add_argument("-i", "--input", help="Input directory containing FASTA files on which the BLAST was run. All FASTA files will be processed.", type=str)
    input_mutex_group.add_argument("--custom-prot2org", help="Custom protein-to-organism mapping file (CSV with columns: protein ID, organism, but no header). If not provided, mapping will be determined from FASTA headers or file names based on --orgbyfile flag.", type=str)

    input_output_group = run_parser.add_argument_group("PSI-BLAST Input and Output", description="Required PSI-BLAST input and output paths")
    input_output_group.add_argument("-u", "--unfiltered_blast", required=True, help="Directory containing unfiltered PSI-BLAST .tab files", type=str)
    input_output_group.add_argument("-f", "--filtered_blast", required=True, help="Directory containing filtered PSI-BLAST .tab files", type=str)
    input_output_group.add_argument("-o", "--outdir", required=True, help="Output directory for all results", type=str)
    
    algorithm_group = run_parser.add_argument_group("Algorithm Parameters", description="Parameters for the triangle clustering algorithm")
    algorithm_group.add_argument("-e", "--evalue_threshold", type=float, default=10.0, help="E-value threshold for BLAST filtering (default: 10.0)")
    algorithm_group.add_argument("-t", "--overlap_threshold", type=overlap_range, default=0.75, help="Overlap threshold (default: 0.75)")
    algorithm_group.add_argument("--mode", choices=["multi","single"], default="multi", help="Allow a protein to be in more than one COG (multi) or keep only occurrence in biggest COG (single)")
    algorithm_group.add_argument("--job", default=None, help="Optional job definition file specifying query and outgroup organisms (format: query_org,outgroup_org).\nIf omitted, each organism is compared against all other organisms (default all-vs-all mode).")
    
    processing_group = run_parser.add_argument_group("Processing Options", description="Options for processing BLAST files and output")
    processing_group.add_argument("-s", "--start_number", type=int, default=1, help="Starting number (default: 1)")
    processing_group.add_argument("-n", "--cog_name", type=str, default="cog")
    processing_group.add_argument("--filter_size", "-fs", type=int, default=1, help="Minimum number of unique organisms required for a cluster to be included in the output. If not set, no filtering is applied.")
    processing_group.add_argument("--format", choices=["legacy","base","minimal","multi"],default="base",nargs="+", help="Output formats:\n" "base = protein name, organism, protein ID, length, start, end, cluster ID;\n""legacy = protein name, organism, protein name, length, start, end, cluster ID;\n" "minimal = protein name, organism, cluster ID;\n" "multi = Collapsed form with only one line per protein and \"multi\" keyword for proteins in more than one cluster")
    processing_group.add_argument("--readblast-output", action="store_true", help="Whether to save intermediate files from readblast step (hits.csv, query2subject.tsv, self.csv, prot2org.csv) in a subdirectory of the output directory")
    processing_group.add_argument("--orgbyfile", action="store_true", help="Whether to determine organism from file name instead of FASTA header.")
    processing_group.add_argument("--suffix", help="File suffix to filter BLAST files", default=".tab")

    run_parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging for debugging purposes, output is saved in \"debug.log\"")

    blast_column_group = run_parser.add_argument_group("BLAST Column Mapping", description="Specify column indices for relevant fields in BLAST files (0-based, default assumes standard BLAST tabular output format 6)")
    blast_column_group.add_argument("-q", "--query_id", help="Column index for query IDs in BLAST files", type=int, default=0, choices=range(12))
    blast_column_group.add_argument("--subject_id", help="Column index for subject IDs in BLAST files", type=int, default=1, choices=range(12))
    blast_column_group.add_argument("--qStart", help="Column index for query start positions in BLAST files", type=int, default=6, choices=range(12))
    blast_column_group.add_argument("--qEnd", help="Column index for query end positions in BLAST files", type=int, default=7, choices=range(12))
    blast_column_group.add_argument("--sStart", help="Column index for subject start positions in BLAST files", type=int, default=8, choices=range(12))
    blast_column_group.add_argument("--sEnd", help="Column index for subject end positions in BLAST files", type=int, default=9, choices=range(12))
    blast_column_group.add_argument("--evalue_id", help="Column index for e-value in BLAST files", type=int, default=10, choices=range(12))
    blast_column_group.add_argument("--score_id", help="Column index for score in BLAST files", type=int, default=11, choices=range(12))

    def run_pipeline(args):
        
        # Clear log file
        with open("debug.log", "w") as f:
            pass
        # Get level
        level = logging.DEBUG if args.verbose else logging.INFO
        setup_logging(level=level)

        # Convert outdir to Path object
        outdir = Path(args.outdir)

        # Warn if directory already exists
        if outdir.exists():
            logger.warning("Warning: Output directory already exists: %s. Files may be overwritten.\n", outdir)

        # Create output directory only after validation
        outdir.mkdir(parents=True, exist_ok=True)

        logger.info("Starting Triangle Clustering Pipeline")
        logger.info("======================================")
        logger.info("Output directory: %s\n", outdir)

        # Define paths for intermediate files
        lse_path = outdir / "lse.csv"

        # ----- Step 1: process input -----
        logger.info("[1/3] Processing input...")
        try:
            # import functions
            from COGprocessinput.map_prot2org import map_prot2org, custom_prot2org_mapping
            from COGprocessinput.readblast import process_blast

            if args.custom_prot2org:
                logger.info("Using custom protein-to-organism mapping from file: %s", args.custom_prot2org)
                prot2org_df = custom_prot2org_mapping(args.custom_prot2org)
                if args.orgbyfile:
                    logger.warning("Warning: --orgbyfile flag is set but will be ignored since a custom protein-to-organism mapping file is provided.")
            else:
                # run prot2org to create mapping dataframe and dictionary
                prot2org_df = map_prot2org(args.input, args.orgbyfile)
            
            entity_dict = dict(zip(prot2org_df.iloc[:, 0], prot2org_df.iloc[:, 1]))

            # Run readblast to process BLAST files and create hits dataframe, query2subject mapping, and self hits dataframe
            (processed_unfiltered_blast, processed_filtered_blast, self_hits) = process_blast(
                entity_dict,
                unfiltered_directory=args.unfiltered_blast,
                filtered_directory=args.filtered_blast,
                query_id=args.query_id,
                qStart_id=args.qStart,
                qEnd_id=args.qEnd,
                subject_id=args.subject_id,
                sStart_id=args.sStart,
                sEnd_id=args.sEnd,
                evalue_id=args.evalue_id,
                evalue_threshold=args.evalue_threshold,
                score_id=args.score_id,
                suffix=args.suffix
            )

            # store intermediate files if requested
            if args.readblast_output:
                blastconv_dir = outdir / "BLASTconv"
                blastconv_dir.mkdir(parents=True, exist_ok=True)

                logger.info("Saving readblast output files to %s", blastconv_dir)
                
                prot2org_df.to_csv(blastconv_dir / "prot2org.csv", index=False, header=["prot_name", "protein", "organism"])

                processed_unfiltered_blast.to_csv(blastconv_dir / "hits.csv", index=False)
                
                with open(blastconv_dir / "query2subject.tsv", "w") as f:
                    for query, targets in processed_filtered_blast.items():
                        f.write(f"{query}\t" + ",".join(map(str, targets)) + "\n")
                
                self_hits.to_csv(blastconv_dir / "self_hits.csv", index=False)
            
            logger.info("Process input completed successfully -> %s\n", outdir)
        except Exception as e:
            logger.error("Error in readblast: %s", e)
            return

        # ----- Step 2: lse -----
        logger.info("[2/3] Running lse...")
        try:
            from COGlse import lse
            presence_map = processed_filtered_blast   # dict[int, list[int]]
            
            df_lse = lse.run_lse_in_memory(
                hits=processed_unfiltered_blast,
                prot2org=prot2org_df,
                presence_map=presence_map,
                output_path=str(lse_path),
                job_path=args.job,
            )
            logger.info("LSE completed successfully -> %s (n_groups=%d)\n", lse_path, len(df_lse))
        except Exception as e:
            logger.error("Error in lse: %s", e)
            return

        logger.debug("LSE DataFrame shape: %s", df_lse.shape)
        


        # ----- Step 3: triangles -----
        logger.info("[3/3] Running triangles...")
        try:
            from COGtriangle.python.cogtriangles import cogtriangles

            cogtriangles(args.evalue_threshold, args.overlap_threshold,processed_unfiltered_blast, prot2org_df, self_hits,
                         df_lse, processed_filtered_blast, args.cog_name, args.start_number, args.outdir, args.mode, args.format, args.filter_size, level)
        except Exception as e:
            logger.error("Error in triangles: %s", e)
            return
    
        logger.info("Pipeline completed successfully.")
        logger.info("All results are stored in: %s\n", outdir)

    run_parser.set_defaults(func=run_pipeline)

    # read arguments
    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
