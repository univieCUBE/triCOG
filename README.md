# triCOG

## Introduction

### What is triCOG?

It is a re-implementation of the triangle-based clustering method for orthologous genes and proteins originally introduced as COGsoft by David Kristensen _et al._ ([2010](https://doi.org/10.1093/bioinformatics/btq229)). triCOG preserves the core methodology while improving usability, modularity, and integration within contemporary software environments.

Please note that triCOG does not reimplement the original COGcognitor functionality, which was designed to incorporate new sequences into pre-existing COG databases. Instead, this project focuses exclusively on the _de novo_ construction of clusters.

### Algorithmic Explanation

The clustering algorithm used in this tool is completely based on the one invented by David Kristensen _et al._ ([2010](https://doi.org/10.1093/bioinformatics/btq229)). It is a deterministic algorithm, originally called EdgeSearch Algorithm.

In the clustering step (`COGtriangles`), the input data is first initialized and then the proteins of each lineage-specific expansion (LSE) are scanned for valid bidirectional hits. Afterwards, the actual clusters are formed: the algorithm picks an unprocessed edge to serve as a seed and then looks for a third vertex that completes a triangle. These new edges are added to a queue, and by repeating this process, cohesive clusters are created from overlapping triangles.

If only two vertices can be combined because no third partner completes a triangle, these are reported as TWOGs (Two-Gene clusters). If no symmetric hits can be found at all, the proteins are returned as singletons.

### Comparison with Similar Tools

Orthology inference is generally split into two methodologies:

- **Tree-based methods (eg. OrthoFinder):** Rely on gene/species tree reconciliation. It distinguishes orthologs from paralogs by annotating duplication vs. speciation events.
- **Graph-based methods (eg. OrthoMCL, InParanoid, COGsoft, or triCOG):** Identify bidirectional best hits (BBHs) between proteins, representing them as edges between the protein nodes. Clusters are then assembled using algorithms such as Markov Clustering (MCL) or symmetric best-hit triangulation (Li _et al._, [2003](https://doi.org/10.1101/gr.1224503); Kristensen _et al._, [2010](https://doi.org/10.1093/bioinformatics/btq229); Emms and Kelly, [2019](https://doi.org/10.1186/s13059-019-1832-y); Persson and Sonnhammer, [2023](https://doi.org/10.1016/j.jmb.2023.168001)).

### When to use triCOG

triCOG is best suited for researchers prioritizing reproducibility and conservative clustering in microbial datasets.

- **Target Data:** Ideal for bacterial, archaeal, and viral proteomes.
- **Reproducibility:** Unlike OrthoMCL (which uses a non-deterministic MCL step), triCOG is fully deterministic. Given the same input, it will always produce the exact same orthologous groups.
- **Precision:** Choose triCOG when you require a conservative estimate (smaller groups with fewer paralogs) and do not need explicit phylogenetic trees (Kristensen _et al._, [2010](https://doi.org/10.1093/bioinformatics/btq229)).

### Limitations

While triCOG provides high-precision clustering, there are scenarios where other tools might be a better fit:

- **Computational Complexity:** triCOG has a worst-case time complexity of $O(n^3 \log n)$ relative to the number of proteomes ($n$). While this is a significant improvement over the original COG algorithm ($O(n^6)$), it can still become a bottleneck for exceptionally large datasets.
- **High Initial Search Cost:** The pipeline relies on pre-computed PSI-BLAST results. This is significantly more computationally expensive than the DIAMOND or BLASTp searches used by other modern orthology tools, especially for large eukaryotic proteomes. However, PSI-BLAST is not required, BLASTp works as well.
- **Eukaryotic Constraints:** For eukaryotic datasets with extensive paralogy or incomplete sequences, tools like OrthoMCL or OrthoFinder are generally preferred. They are specifically benchmarked to handle the complex duplication events found in higher organisms.
- **No Phylogenetic Inference:** triCOG is a clustering tool, not a tree-builder. It does not infer gene or species trees. If your research requires explicit phylogenetic context or reconciliation, consider using SYNERGY or OrthoFinder (Li _et al._, [2003](https://doi.org/10.1101/gr.1224503); Wapinski _et al._, [2007](https://doi.org/10.1093/bioinformatics/btm193); Emms and Kelly, [2019](https://doi.org/10.1186/s13059-019-1832-y)).
- **Memory-trade-offs:** Preliminary benchmarks indicate that triCOG notably reduces execution time compared to COGsoft. However, this speedup comes at the cost of a higher memory footprint. While optimised for modern systems, further memory management refinements are required to ensure scalability for massive datasets (see [Future Considerations](DEVELOPMENT.md#performance-optimizations)).

## Getting started

### Prerequisites

To get started, ensure your environment meets the following requirements:

- **OS:** UNIX-based system (Linux/macOS)
- **Package Manager:** Conda
- **Compiler:** `C++17` compatible compiler (e.g., GCC 7+, Clang 5+)

A full list of dependencies is available in the dependency section of the [Developer Documentation](DEVELOPMENT.md#dependencies).

### Installation

The software can be set up by using the `environment_triCOG.yml` file to create a suitable conda environment:

```bash
conda env create -f environment_triCOG.yml  # Create the conda environment called 'triCOG_env'
conda activate triCOG_env                   # Activate the conda environment.
```

After creating the conda environment, the the `C++` module has to be compiled. This can be done with the `setup_compilation.sh` bash script.

```bash
bash setup_compilation.sh                   # Compiling the program
```

### Verification

The following lines can be run to correctly test if the compilation worked:

```bash
python triCOG.py -h                         # Should print the general help message
python triCOG.py run -h                     # Should print the specific run pipeline help message
```

## Tutorial

### Preparing input

#### BLAST results

Provide two separate folders containing your BLAST results, PSI-BLAST ([Altschul and Koonin, 1998](<https://doi.org/10.1016/S0968-0004(98)01298-5>)) is recommended. The folders can contain multiple BLAST output files, but all of them need to have the same file suffix (default: `.tab`).

- `--unfiltered_blast`: Results generated without composition-based statistics (Schäffer _et al._, [2001](https://doi.org/10.1093/nar/29.14.2994)) and no SEG filtering (Wootton and Federhen, [1996](<https://doi.org/10.1016/S0076-6879(96)66035-2>)).
- `--filtered_blast`: Results generated with composition-based score adjustment (Schäffer _et al._, [2001](https://doi.org/10.1093/nar/29.14.2994)) and SEG filtering (Wootton and Federhen, [1996](<https://doi.org/10.1016/S0076-6879(96)66035-2>)).

#### Organism Input

One of the following ways to map proteins to organisms has to be chosen:

- `--input`: A directory containing the original FASTA files (`.fa`, `.fasta`, `.fna`, `.faa`). The directory content can be single or multi files. Currently, NCBI [three-letter tag](https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5) for the protein sequence identifier (SeqID) are [not supported](DEVELOPMENT.md#feature-extensions).
- `--custom-prot2org`: A custom prot2org CSV file containing `<protein_name>,<organism_name>` pairs that should be consider for the clustering. The `protein_name` needs to match the one in the PSI-BLAST result. `protein_name` is generally the SeqID.

#### LSE job file

The LSE module handles query-outgroup relationships. You can either let the system automate this or provide a manual configuration.

| Mode      | Flag         | Behavior                                                                                     |
| --------- | ------------ | -------------------------------------------------------------------------------------------- |
| Automatic | (None)       | All-vs-all: Each organism is treated as a query, with all others acting as outgroups.        |
| Manual    | `--job-file` | Custom: Define specific biological groups or restrict the analysis to a subset of organisms. |

**Job File Format** <br>
Each line must define at least one query-outgroup pair: <br>`<query_organism>,<outgroup_organism>`

### Running the analysis

The complete pipeline can be run using the following command:

```bash
python triCOG.py run [options] -i INPUT_FASTA_FOLDER -u UNFILTERED_BLAST -f FILTERED_BLAST -o OUTDIR
```

For the specific options, see the [required argument section](#required-arguments).

#### Modular execution

The three separate parts of the program, COGprocessinput, COGlse and COGtriangle can be run separately via subcommands.
The complete pipeline is like running the three subcommands one after the other with the output from the previous subcommand as the input of the next. The modules can be run separately like this:

```bash
python triCOG.py processinput [options] -i INPUT_FASTA_FOLDER -u UNFILTERED_BLAST -f FILTERED_BLAST -o OUTDIR

python triCOG.py lse [options] --hits OUTDIR/hits.csv --prot2org OUTDIR/prot2org.csv --q2s OUTDIR/query2subject.tsv

python triCOG.py triangles [options] --hits OUTDIR/hits.csv --prot2org OUTDIR/prot2org.csv --prot2len OUTDIR/self_hits.csv --q2s OUTDIR/query2subject.tsv --lse lse.csv --output OUTDIR
```

Specific arguments for the subcommands can be found in the [arguments for modular execution](#arguments-for-modular-execution) section.

### Output

The output consists of the following files:

| File                  | Content                                                                                             | Structure                                                                           |
| --------------------- | --------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `{format}_output.csv` | The output, where every protein contained in the input fasta files is assigned a cluster (or none). | Format depends on `--format` flag (specified in [output formats](#output-formats)). |
| `all_edges.csv`       | All edges used to create the clusters.                                                              | Each line is a pair of 2 protein SeqIDs and their hash integer which form the edge. |
| `cog_edges.csv`       | All edges assigned to the specific cluster they form.                                               | Each line is an edge, with added cluster ID and organism name.                      |
| `lse.csv`             | All Lineage-specific expansions.                                                                    | Each line is one LSE group.                                                         |

The tool offers several post-processing steps to refine the clustering result. Besides the choice of different output formats, the functionality of the original `COGtriangles.reformat.pl` was implemented.

#### Clustering modes

Users can choose between two primary output strategies:

- **Multi Cluster Mode (default):** This is what the algorithm generates. Proteins are allowed to be members of multiple clusters simultaneously.
- **Single Cluster Mode:** Ensures each protein appears only once. If a protein is found in multiple clusters, the algorithm keeps only the occurrence within the largest cluster.

#### Output Formats

The tool provides four output formats. In all formats, each row represents a protein and its associated cluster. The output structure follows the original COGsoft convention: the `start` coordinate is hardcoded to `1`, and the `end` coordinate corresponds to the full protein length. According to the COGsoft documentation, this design allows for advanced use cases in which proteins may be subdivided into domains with `start` and `end` positions defined relative to the full-length sequence. However, in the current implementation, proteins are always treated as full-length sequences (see [Future Considertations](DEVELOPMENT.md#feature-extensions)).

Available Formats

- `base`: default format <br>
  protein SeqID, organism name, protein hash integer, alignment length, 1, alignment length, cluster ID
- `legacy`: identical to COGsoft implementation <br>
  protein SeqID, organism name, protein SeqID, alignment length, 1, alignment length, cluster ID
- `minimal`: no display of BLAST information
  protein SeqID, organism name, cluster ID
- `multi`: only available in `multi mode` <br>
  collapses the content to only one line per protein analogous to the `-b` mode of the COGsoft tool

#### Filtering size

The `--filter-size` argument filters the output based on the number of organisms in a cluster. It can be set to any integer. The most relevant are:

- `--filter-size 1`: No filtering is performed (default).
- `--filter-size 2`: Singletons (proteins without a cluster) are removed.
- `--filter-size 3`: Only triangle based clusters remain.

For massive datasets, where only triangle based clusters are needed, the processing time can be reduced by disabling the creation of Singletons and TWOGs at compiler level. If this optimisation is enabled, the program should still be called with `--filter-size 3` to ensure correct results for the `single COG mode` and `multi` output.

To recompile, use this command at the root directory:

```bash
# Recompile with optimisation flag
cmake -S . -B build -DNO_TWOG_SINGLETON=ON
cmake --build build -j
```

## Examples

### Required arguments only

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/
```

### Custom thresholds

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/ -e 0.01 -t 0.5
```

### Single COG mode

Proteins occur only in the biggest cluster (see [clustering modes](#clustering-modes)).

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/ --mode single
```

### Organism assignment by file name

Organism is not chosen per fasta definition line modifier, but by the file name the sequence is in.

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/ --orgbyfile
```

### Custom prot2org file

```bash
python triCOG.py run --custom-prot2org prot2org_path  -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/
```

### Custom job file for the LSE

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/ --job job_filepath
```

### Filter for triangle clusters

Keep only clusters with at least 3 different organisms in output.

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/ --filter-size 3
```

### Legacy output

Output format identical to COGsoft implementation by Kristensen _et al._ ([2010](https://doi.org/10.1093/bioinformatics/btq229)).

```bash
python triCOG.py run -i fasta_folder/ -u unfiltered_blast/ -f filtered_blast/ -o COGoutput/ --format legacy
```

## Command-Line references

### Required Arguments

There are five required arguments, listed in the table below. `--unfiltered_blast`, `--filtered_blast` and `--outdir` are always required, while `--input` and `--custom-prot2org` are mutually exclusive. `--custom-prot2org` is recommended when the FASTA definition line does not contain the organism modifier and the FASTA files are not separated by organism.

| Argument                 | Description                                                                                                                                                                                                   | Type           |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------- |
| `-i , --input`           | Input directory containing FASTA files on which the BLAST was run. All FASTA files will be processed.                                                                                                         | Directory path |
| `--custom-prot2org`      | Custom protein-to-organism mapping file (CSV with columns: protein SeqID, organism, but no header). If not provided, mapping will be determined from FASTA headers or file names based on `--orgbyfile` flag. | File Path      |
| `-u, --unfiltered_blast` | Directory containing unfiltered PSI-BLAST output files. All files will be processed.                                                                                                                          | Directory path |
| `-f, --filtered_blast`   | Directory containing filtered PSI-BLAST output files. All files will be processed.                                                                                                                            | Directory path |
| `-o, --outdir`           | Output directory for all output files. Existing files can be overwritten.                                                                                                                                     | Directory path |

### Algorithm Parameters

| Argument                  | Description                                                                                                                                    | Default | Type           |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | ------- | -------------- |
| `-e, --evalue_threshold`  | E-value threshold for BLAST filtering                                                                                                          | 10.0    | Float          |
| `-t, --overlap_threshold` | Overlap threshold (range 0-1)                                                                                                                  | 0.75    | Float          |
| `--mode`                  | Allow a protein to be in multiple COGs (`multi`) or keep only occurrence in biggest COG (`single`)                                             | `multi` | String         |
| `--job`                   | Optional job file specifying query and outgroup organisms (format: query_org,outgroup_org). If omitted, default all-vs-all outgroups are used. | `-`     | Directory Path |

### Processing Options

| Argument             | Description                                                                                                                                                                | Default | Type    |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------- | ------- |
| `-s, --start_number` | Start number for COG naming                                                                                                                                                | `1`     | Integer |
| `-n, --cog_name`     | Prefix for COG naming                                                                                                                                                      | `cog`   | String  |
| `--filter_size`      | Filter output for clusters with at least this many organisms per cluster.                                                                                                  | `1`     | Integer |
| `--format`           | Output format (`base`, `legacy`, `minimal` or `multi`)                                                                                                                     | `base`  | String  |
| `--readblast-output` | Whether to save intermediate files from processinput step (hits.csv, query2subject.tsv, self.csv, prot2org.csv) in a subdirectory of the output directory called BLASTconv | `-`     | Flag    |
| `--orgbyfile`        | Whether to determine organism from file name instead of FASTA header.                                                                                                      | `-`     | Flag    |
| `--suffix`           | File suffix of BLAST files. Only files with this suffix are processed.                                                                                                     | `.tab`  | String  |

### BLAST input column mapping

The expected output format is the BLAST output format 6 or 7. The default values match those output formats. For other formats, the columns can be specified as integers (0-based) with the following arguments:

| Argument         | Description                                             | Default | Valid range |
| ---------------- | ------------------------------------------------------- | ------- | ----------- |
| `-q, --query_id` | Column index for query SeqIDs in BLAST files.           | `0`     | 0-11        |
| `--subject_id`   | Column index for subject SeqIDs in BLAST files.         | `1`     | 0-11        |
| `--qStart`       | Column index for query start positions in BLAST files   | `6`     | 0-11        |
| `--qEnd`         | Column index for query end positions in BLAST files     | `7`     | 0-11        |
| `--sStart`       | Column index for subject start positions in BLAST files | `8`     | 0-11        |
| `--sEnd`         | Column index for subject end positions in BLAST files   | `9`     | 0-11        |
| `--evalue_id`    | Column index for e-value in BLAST files                 | `10`    | 0-11        |
| `--score_id`     | Column index for bit score in BLAST files               | `11`    | 0-11        |

All columns mentioned above are mandatory for the program to work.

### Other Arguments

| Argument        | Description                                                                   | Default | Type |
| --------------- | ----------------------------------------------------------------------------- | ------- | ---- |
| `-v, --verbose` | Enable verbose logging for debugging purposes, output is saved in "debug.log" | `-`     | Flag |
| `-h, --help`    | Show help message and exit.                                                   | `-`     | Flag |

## Arguments for modular Execution

While `triCOG` offers an end-to-end pipeline, every stage is available as a standalone module. This decoupled approach is ideal for debugging specific steps or for scenarios requiring high-precision manual control. In standalone mode, each module operates independently by reading and writing directly to your filesystem.

### processinput Subcommand

The `processinput` subcommand takes the [same input files](#required-arguments) as the complete pipeline. Additional options are `-e`, `--orgbyfile`, `--custom-prot2org`, and `--suffix` from the [processing option table](#processing-options) and all [BLAST column mappings](#blast-input-column-mapping). It is called by the command:

```bash
python triCOG.py processinput [options] -i INPUT_FASTA_FOLDER -u UNFILTERED_BLAST -f FILTERED_BLAST -o OUTDIR
```

### lse subcommand

The `lse` subcommand takes additional input arguments compared to the full run. All the needed data structures have to be specified as:

| Argument     | Description                                                                                           | Default   | Type           |
| ------------ | ----------------------------------------------------------------------------------------------------- | --------- | -------------- |
| `--hits`     | The generated hits table from COGprocessinput `hits.csv`                                              | `-`       | Directory Path |
| `--prot2org` | Hashing and organism file (protein SeqID, hash integer, organism name)                                | `-`       | Directory Path |
| `--q2s`      | Filter TSV file to limit possible edges (query, subjects), e.g. COGprocessinput's `query2subject.tsv` | `-`       | Directory Path |
| `--output`   | Output file                                                                                           | `lse.csv` | String         |

Moreover, the [`--job` algorithm parameter](#algorithm-parameters) and `--verbose` can be used as well.

### triangles Subcommand

The `triangles` subcommand takes additional input arguments compared to the full run. All the needed data structures have to be specified as:

| Argument     | Description                                                            | Default | Type           |
| ------------ | ---------------------------------------------------------------------- | ------- | -------------- |
| `--hits`     | The generated hits table from COGprocessinput `hits.csv`               | `-`     | Directory path |
| `--prot2org` | Hashing and organism file (protein SeqID, hash integer, organism name) | `-`     | Directory path |
| `--prot2len` | Protein-to-length mapping file (name, length)                          | `-`     | Directory path |
| `--lse`      | Lineage-specific expansion file (one line per LSE)                     | `-`     | Directory path |
| `--q2s`      | Filter file to limit possible edges (query, subjects)                  | `-`     | Directory path |

Moreover, the [algorithm parameters](#algorithm-parameters) and `--start_number`, `--cog_name`, `--filter_size`, `--format`, `--outdir` and `--verbose` can be used to specify the clustering algorithm and output format.

## References

Altschul, S. F., & Koonin, E. V. (1998). Iterated profile searches with PSI-BLAST—a tool for discovery in protein databases. Trends in biochemical sciences, 23(11), 444-447.

Emms, D. M., & Kelly, S. (2019). OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome biology, 20(1), 238.

Kristensen, D. M., Kannan, L., Coleman, M. K., Wolf, Y. I., Sorokin, A., Koonin, E. V., & Mushegian, A. (2010). A low-polynomial algorithm for assembling clusters of orthologous groups from intergenomic symmetric best matches. Bioinformatics, 26(12), 1481-1487.

Li, L., Stoeckert, C. J., & Roos, D. S. (2003). OrthoMCL: identification of ortholog groups for eukaryotic genomes. Genome research, 13(9), 2178-2189.

Persson, E., & Sonnhammer, E. L. (2023). InParanoiDB 9: ortholog groups for protein domains and full-length proteins. Journal of molecular biology, 435(14), 168001

Schäffer, A. A., Aravind, L., Madden, T. L., Shavirin, S., Spouge, J. L., Wolf, Y. I., ... & Altschul, S. F. (2001). Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics and other refinements. Nucleic acids research, 29(14), 2994-3005.

Wapinski, I., Pfeffer, A., Friedman, N., & Regev, A. (2007). Automatic genome-wide reconstruction of phylogenetic gene trees. Bioinformatics, 23(13), i549-i558.

Wootton, J. C., & Federhen, S. (1996). [33] Analysis of compositionally biased regions in sequence databases. In Methods in enzymology (Vol. 266, pp. 554-571). Academic Press.
