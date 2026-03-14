# triCOG

## Developer-level Documentation

This project is a re-implementation of the triangle-based clustering method for orthologous genes and proteins originally introduced as COGsoft by David Kristensen _et al._ ([2010](https://doi.org/10.1093/bioinformatics/btq229)). triCOG preserves the core methodology while improving usability, modularity, and integration within contemporary software environments.

The software is organized into three main components:

- COGprocessinput – Preprocessing and preparation of input data
- COGlse – Computation of lineage-specific expansion (LSE)
- COGtriangles – Execution of the triangle-based clustering algorithm

The core application logic is implemented in `Python`, providing flexibility and ease of integration. The computationally intensive clustering algorithm has been preserved and optimized in `C++` and is exposed to `Python` through a `pybind11` interface to ensure high performance while maintaining a clean `Python` API.
A tutorial on how to use the software can be found in the [User Documentation](README.md).

## Architectural Overview

The command-line interface is defined in `triCOG.py` in the project’s root directory. This file serves as the main entry point and coordinates program execution based on the provided options.

### COGprocessinput

The purpose of COGprocessinput is, as the name suggests, to process the input data to fit the input data structures required by [COGlse](#coglse) and [COGtriangles](#cogtriangles). The [`map_prot2org.py` module](COGprocessinput/map_prot2org.py) consists of the `map_prot2org(input_directory, OrgByFile)` function, which reads all protein sequence identifiers (SeqID) from the fasta files specified in the `-i`/`--input` command line argument, and assigns SeqID a hashing value and associates it with its organism. This can either be done via the file name where the protein sequence is in (when the `--orgbyfile` flag is specified) or by reading the modifiers in the fasta definition line to find the species name there, which is explained in more detail in the [Implementation Details](#automatization-of-prot2org-creation). Alternatively, it can also read in a custom `prot2org` file. To avoid confusion between the protein sequence ID and the hash ID of each protein, we use the names `prot_name` or `protein_name` for the protein SeqID, and the name `protein` for the hash integer. The resulting data structure is a pandas DataFrame which consists of three columns: protein name, hash integer, organism name.

Additionally, the [`readblast.py` module](COGprocessinput/readblast.py) contains the `process_blast()` function, which iterates through both the directory containing the unfiltered PSI-BLAST results and filtered PSI-BLAST results. In this context, unfiltered means that no composition-based statistics and no SEG filtering were used, while for the filtered PSI-BLAST run composition-based score adjustment and SEG filtering was used (Wootton and Federhen, [1996](<https://doi.org/10.1016/S0076-6879(96)66035-2>); Schäffer _et al._, [2001](https://doi.org/10.1093/nar/29.14.2994)). Each BLAST hit in the directories is read and checked for completeness and whether it passes the e-value threshold. The default values for the parsed columns are the output formats 6 or 7 of BLAST. Hits from the unfiltered BLAST directory are stored in a pandas DataFrame called `unfiltered_blast_results`, while for hits from the filtered BLAST directory only the query and subjects are stored in a dictionary called `filtered_blast_results`, where the keys are queries and the values are lists of all subjects the key had a hit with. Based on `unfiltered_blast_results`, a third data structure is created, which is either a self hit for each query, or if the query does not have a self hit, then it takes the hit with the best bit score for that query. These self hits are saved in a pandas DataFrame called `self_hits`. These three data structures are then used by both COGlse and COGtriangles, together with the `prot2org` DataFrame.

### COGlse

The COGlse component computes Linked Sequence Elements (LSEs) from BLAST hits and a protein-to-organism mapping. The implementation is located in the [`lse.py` module](COGlse/lse.py). In the `Python` pipeline, the LSE step uses the unfiltered BLAST hits as input (`processed_unfiltered_blast`) and applies strict presence constraints from the filtered BLAST representation (`processed_filtered_blast`).

Processing is performed organism-wise. For each query organism, all its proteins are considered as nodes. Hits are first annotated with organism labels and then sorted by query, then `e-value` (ascending) and `score` (descending). For each query protein, hits are scanned in this order. As soon as a hit to an outgroup organism (specified in the `job file`) is found, scanning stops for this query protein. Only intra-genome hits that appear before this outgroup hit are kept.

From these directed intra-genome hits, reciprocal connections are extracted. An undirected edge between two proteins is created only if both directions (a -> b and b -> a) exist. These edges are then filtered using the `query2subject` map. An edge is kept only if both proteins appear in each other’s presence list (`query2subject`).
For the current organism, an undirected graph is built using all its proteins as nodes and the filtered reciprocal edges as connections. Connected components of this graph define the final LSE groups. Proteins without reciprocal connections form single member groups.

All groups across organisms are collected into a single DataFrame (`lse_df`). This DataFrame is returned to the pipeline for downstream processing (triangle clustering) and can optionally be written to disk (single-column, no header) as `lse.csv` file.

### COGtriangles

The actual clustering of the data is performed in the COGtriangles component. The [`cogtriangles.py` module](COGtriangle/python/cogtriangles.py) (located in the `python/` subdirectory) coordinates data flow between `Python` and the `C++` backend. Input data structures are converted into `Python` dictionaries to enable automatic type conversion via `pybind11` in the [`data_converter.py`](COGtriangle/python/data_converter.py) file. In accordance with the original implementation, the data is then initialized as a `HitSet` class instance in `C++`. After that the `C++` `insert()` function is called for every LSE in the given dataset. In the function, the best hit is found for each protein by checking the hits table. Here, the `query2subject` filter is used to decide whether a protein is included in the search for best hit or not. The best hits for the proteins are then collapsed on a LSE level for further processing. If one target has more than half of the hits of the current LSE this LSE is chosen as best hit, else if it is exactly half, the further proceeding is handled separately in the conflict map. If no target reaches 50 % of the hits, no best hit is found.

After all proteins have been processed, [`cogtriangles.py`](COGtriangle/python/cogtriangles.py) invokes the `makeCOGs()` function to compute the final clusters from the constructed edge set. Here, all bidirectional edges are saved in `all_edges.csv` and the triangle algorithm is seeded with unused edges. The seed is then extended with a breadth-first-traversal by trying to find edges that complete a triangle of the existing cluster. All edges that are part of clusters are saved in the `cog_edges.csv`. Two organisms clusters and singletons are handled separately and added to the output.

The resulting clusters are returned to `Python` for post-processing. Protein sequence identifiers are restored, and additional processing steps are applied in [`process_output.py`](COGtriangle/python/process_output.py). This module handles the distinction between `single-COG` and `multi-COG` modes, determines the output format, and applies optional size-based filtering (more details in [README](README.md#output-formats) and [Implementation Details](#output-processing)). Final results are written to the user-specified output directory.

## Implementation Details

### `C++` integration

Before using the CLI tool, the `C++` clustering module must be compiled. We provide a standalone shell script ([`setup_compilation.sh`](setup_compilation.sh)) for this purpose to ensure a flexible build process across different environments. For the compilation a compiler with `C++17` support (GCC 9+ or Clang 5+) is needed.

### CLI

The command-line interface leverages `Python`’s built-in argparse module to parse and validate user-provided arguments. For improved debugging and flexibility, each of the three components can also be executed independently as standalone commands, using the syntax described in the [User Documentation](README.md).

### Automatization of prot2org creation

The original COGsoft implementation (Kristensen _et al._, [2010](https://doi.org/10.1093/bioinformatics/btq229)) required a mapping of each protein SeqID (`prot_name`) to an organism name (here called `prot2org`), which the user was supposed to create beforehand. To automatize this step, the [`map_prot2org.py` module](COGprocessinput/map_prot2org.py) can now iterate through all fasta files in a directory, ideally the same fasta files used to run the PSI-BLAST search and do this mapping.

There are two ways for this mapping. The first way is that each protein SeqID is mapped to the file it is in, indicated by the `--orgbyfile` flag. The other way is to extract the species name of each protein SeqID from its fasta description line. For that, the assumption is that the fasta description line follows the [NCBI guideline for modifiers](https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/) in a fasta header, i.e. the species is either clarified via the `organism` modifier or be the rightmost modifier in the definition line without an explicit modifier.

### Makehash Replacement

The first step of the original COGsoft implementation was to create a hashing for all proteins in the `prot2org` mapping, which was done by a subprogram called COGmakehash. This function was taken over by the [`map_prot2org.py`](COGprocessinput/map_prot2org.py) [module](#automatization-of-prot2org-creation), and was therefore not required anymore. The hashing itself is still done, as working with integers is still the faster and more memory-efficient way. Memory usage was further reduced by combining the SeqID-to-hash and SeqID-to-organism mapping into the `prot2org` DataFrame.

### Change `query2subject` data structure

The `query2subject` data structure was changed from a DataFrame in the original COGsoft implementation (Kristensen _et al._, [2010](https://doi.org/10.1093/bioinformatics/btq229)) to a dictionary in the form of `{query: [subject1, subject2, ...]}` (with query and subject as their respective hash integers). This was done to eliminate storing each query multiple times to reduce memory usage, and it also speeds up retrieving and iterating over all subjects of a query with the hash-based lookup of dictionaries, compared to DataFrame column filtering or row-wise iteration.

### LSE job file

In contrast to the original software, if no job file is provided, the module automatically generates a default job configuration (`generate_default_job`, `get_job`). In this configuration, each organism is processed as a query organism while all remaining organisms are treated as outgroups. This allows the LSE module to run without requiring an explicit job file and simplifies usage in standard all-vs-all scenarios. At the same time, a custom job file can be provided to define specific query-outgroup relationships. A job file must at least define one query-outgroup pair per line in the format: `<query_organism>,<outgroup_organism>`

### Edge Filtering (query2subject Logic)

The `query2subject.tsv` file, generated during the COGprocessinput phase, acts as a whitelist for graph edges. The algorithm implements a conditional filtering logic:

- **Empty `query2subject`:** The filter is bypassed and all hits that pass the thresholds are processed.
- **Non-empty `query2subject`:** Only the edges explicitly defined in this file are permitted. All other potential hits are discarded.

### Logging Strategy

The software uses two separate logging systems to handle the `Python`-`C++` boundary. While a single unified logger would be ideal, a decoupled approach was chosen to ensure stability and avoid complex cross-language stream handling and uninitialized loggers.

- **`Python` side:** Uses the native `logging` module. It captures all high-level logic and orchestrates the log level and debug log.
- **`C++` side:** Uses `spdlog` for high-performance logging. It is initialized within the `pybind11` bindings, and the active log level is explicitly passed to the `C++` layer via the `HitSet` to ensure consistent logging behavior across both environments.

For simplified configuration, log output is written to debug.log, which is only generated when debug mode (`-v`) is enabled.

### Output Processing

One of the primary goals of this re-implementation was to keep the same functionality as the underlying tool. Therefore, we decided to integrate `COGtriangles.reformat.pl` logic directly into this tool, eliminating the need for external `Perl` dependencies. Hence, the strict and b mode besides the filtering logic were integrated in this software.

To ensure compatibility with existing pipelines, this tool offers the `legacy output`, which exactly recreates the original output. To include less redundant formats the `base` and `minimal` outputs were added for completion. Moreover, the `-b` mode from the original was implemented as `multi` output.

## Dependencies

The used dependencies can be found in `environment_triCOG.yml`. The following table shows the used version and the suggested minimal needed version of the packages.

| Tool       | Used Version | Suggested Min. Version |
| ---------- | :----------: | :--------------------: |
| `python`   |    3.13.0    |          3.10          |
| `numpy`    |    2.4.0     |          2.0           |
| `pandas`   |    2.3.3     |          2.0           |
| `pybind11` |    3.0.1     |          2.13          |
| `cmake`    |   >= 3.20    |          3.20          |
| `spdlog`   |    1.16.0    |          1.10          |

## Testing

To ensure the reliability of the core `C++` implementation and its `Python` bindings, this project includes a test suite located in the `tests/` directory. The suite is partitioned by component for better maintainability and navigation:

- **`Python` Tests:** Unit and integration tests are implemented using pytest. To run the tests from the project root, use:

```bash
pytest tests/
```

- **`C++` Tests:** Native module tests leverage the `Catch2` framework. The framework is automatically fetched via CMake's FetchContent during the configuration step. After building the project, you can execute the `C++` unit tests directly:

```bash
./build/COGtriangle/unit_tests
```

All [testdata](tests/testdata/) is build upon publicly available proteomes, for which the accession numbers can be found in [`accession_numbers.txt`](tests/testdata/accession_numbers.txt).

## Future Considerations

The current version provides a modernized foundation for clustering, but several areas remain open for future enhancement:

### Performance Optimizations

- **Parallelization:** The `C++` backend currently operates on a single thread. The edge insertion process, in particular, has potential for multi-threading, which would speed up the processing of large-scale datasets.
- **Efficient Data Bridging:** While `pybind11`'s automatic dictionary conversion is robust, transitioning to `NumPy` buffer protocols would enable zero-copy data transfer. This would eliminate the overhead of type conversion and reduce memory consumption.
- **Profiling & Scaling:** Due to current limitations in the testing data, it is difficult to predict exactly how the tool will behave with massive proteomes.
  - **Memory Management:** To make the pipeline more robust, future versions should likely implement chunked data processing. This would avoid loading everything into memory at once, which is currently a potential point of failure for very large datasets.
  - **Future Testing:** Further profiling is needed to identify other unforeseen bottlenecks. It remains to be seen how the current logic scales when hit sets grow by orders of magnitude.

### Feature Extensions

- **Labeling of protein SeqIDs:** The NCBI allows for [three-letter tags](https://ncbi.github.io/cxx-toolkit/pages/ch_demo#ch_demo.T5) for additional information on the SeqID. These are currently not supported, but could be included in the future for a wider use-case range.
- **Output coordinates:** The current implementation follows the original COGsoft convention, where output rows use hardcoded coordinates (`start = 1`, `end = protein length`). While this preserves compatibility, it does not reflect the actual alignment or hit coordinates.
  A future revision could replace these with real alignment-based coordinates (e.g., BLAST hit start/end positions) or optionally support domain-level coordinates. This would make the output more informative and better suited for domain-aware or region-specific analyses.
- **Advanced Output Formats:** To improve integration into other pipelines, output support could be extended beyond `CSV` files to `JSON`, `SQLite` or `Parquet` formats.

## References

Kristensen, D. M., Kannan, L., Coleman, M. K., Wolf, Y. I., Sorokin, A., Koonin, E. V., & Mushegian, A. (2010). A low-polynomial algorithm for assembling clusters of orthologous groups from intergenomic symmetric best matches. Bioinformatics, 26(12), 1481-1487.

Schäffer, A. A., Aravind, L., Madden, T. L., Shavirin, S., Spouge, J. L., Wolf, Y. I., ... & Altschul, S. F. (2001). Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics and other refinements. Nucleic acids research, 29(14), 2994-3005.

Wootton, J. C., & Federhen, S. (1996). [33] Analysis of compositionally biased regions in sequence databases. In Methods in enzymology (Vol. 266, pp. 554-571). Academic Press.
