# lse.py
from __future__ import annotations
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple, Optional, Any
import pandas as pd
from collections import defaultdict
import logging

logger = logging.getLogger("lse")
logging.basicConfig(level=logging.INFO)


# -----------------------------
# 1.  Input / Parsing handling
# -----------------------------

def load_hits(hits_path: str | Path) -> pd.DataFrame:
    """
    Load BLAST hits from disk (CSV format).

    Used by the file-based LSE entry point. Performs minimal validation
    and type normalization required for early-stop sorting.

    Args:
        hits_path: Path to the BLAST hits CSV file
                   (must contain columns: query, subject, score, evalue).

    Returns:
        pd.DataFrame:
            DataFrame with columns ["query", "subject", "score", "evalue"],
            where query/subject are int and score/evalue are float.
    """

    # Allow both string and Path input
    hits_path = Path(hits_path)

    # Read BLAST hits
    df = pd.read_csv(hits_path, sep=",", comment="#")

    # Ensure required columns exist
    required = {"query", "subject", "score" , "evalue"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"hits-file is missing the column: {sorted(missing)}")

    # Ensure that query and subject IDs are integers
    df["query"] = df["query"].astype(int)
    df["subject"] = df["subject"].astype(int)

    # Ensure score and evalue are floats (for correct early-stop sorting)
    df["score"] = df["score"].astype(float)
    df["evalue"] = df["evalue"].astype(float)

    return df


def load_prot2org(prot2org_path: str | Path) -> pd.DataFrame:
    """
    Load protein→organism mapping from disk (CSV format).

    Used by the file-based LSE entry point. Performs minimal validation
    and type normalization required for downstream grouping.

    Args:
        prot2org_path (str | Path): Path to the protein-to-organism CSV file
                                    (must contain columns: protein, organism).

    Returns:
        pd.DataFrame:
            DataFrame with columns ["protein", "organism"],
            where protein is int and organism is a stripped string.
    """

    prot2org_path = Path(prot2org_path)

    # Read mapping table 
    df = pd.read_csv(prot2org_path, sep=",", comment="#",)

    # Ensure required columns exist
    required = {"protein", "organism"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"prot2org-file is missing the column: {sorted(missing)}"
        )

    # ID needs to be int
    df["protein"] = df["protein"].astype(int)

    # Organism needs to be string
    df["organism"] = df["organism"].astype(str).str.strip()

    return df[["protein", "organism"]].copy()


def load_job_file(job_path: str | Path) -> Dict[str, Set[str]]:
    """
    Load LSE job description from disk.

    Used by the file-based LSE entry point to define, for each query organism,
    the set of organisms that should be treated as outgroups.

    Args:
        job_path (str | Path): Path to the job description file.

    Returns:
        Dict[str, Set[str]]: Mapping from query organism to outgroup organisms.
    """

    job_path = Path(job_path)

    # Robust parsing: accept comma or tab separators, ignore comments, no header assumed
    df = pd.read_csv(
        job_path,
        sep=None,              
        engine="python",
        header=None,
        comment="#",
    )

    # Job file must define at least (query_org, outgroup_org)
    if df.shape[1] < 2:
        raise ValueError(
            "job-file must have at least 2 columns: <query_org>,<outgroup_org>"
        )

    df = df.iloc[:, :2].copy()
    df.columns = ["query_org", "outgroup_org"]
    df["query_org"] = df["query_org"].astype(str).str.strip()
    df["outgroup_org"] = df["outgroup_org"].astype(str).str.strip()

    # Build mapping: one query organism can have multiple outgroups
    job: Dict[str, Set[str]] = {}
    for q, og in df.itertuples(index=False):
        if q == "" or og == "":
            continue
        job.setdefault(q, set()).add(og)

    return job


def load_presence_map(path: str | Path) -> Dict[int, List[int]]:
    """
    Load query-to-subject mapping (q2s) from disk.

    In the CLI (triCOG.py) this file is provided via --q2s. Internally it is referred to
    as `presence_map`, since it defines allowed adjacency/presence relations
    between proteins.

    Expected format (TSV, no header):
        query_id<TAB>subject1,subject2,...

    Returns:
        dict[int, list[int]]:
            Mapping query_id -> list of subject_ids.
    """
    path = Path(path)

    q2s_df = pd.read_csv(path, sep="\t", header=None)

    presence = {
        int(k): [int(x) for x in str(v).split(",") if x]
        for k, v in zip(q2s_df.iloc[:, 0], q2s_df.iloc[:, 1])
    }

    return presence



# -----------------------------
# 2.  Job / Setup
# -----------------------------

def build_id_to_org(prot2org_df: pd.DataFrame) -> Dict[int, str]:
    """
    Build protein-ID → organism mapping for the LSE algorithm.

    This converts the prot2org DataFrame into a plain dictionary, which is the
    only representation used by the LSE pipeline from this point on.
    """

    mapping = dict(zip(prot2org_df["protein"], prot2org_df["organism"]))
    return mapping


def list_organisms(id_to_org: Dict[int, str]) -> List[str]:
    """
    List all organism identifiers present in the dataset.

    This defines the set of query organisms for which LSE groups
    will be computed.
    """

    organisms = sorted(set(id_to_org.values()))
    return organisms


def generate_default_job(organisms: List[str]) -> Dict[str, Set[str]]:
    """
    Generate a default LSE job definition.

    A job specifies for which query organism LSE is computed and which
    organisms are treated as outgroups for comparison.

    By default, for each organism all other organisms are used as outgroups.

    Args:
        organisms (List[str]): List of organism identifiers.

    Returns:
        Dict[str, Set[str]]: Mapping from query organism to its set of outgroups.
    """

    job = {}

    # For each organism, define all remaining organisms as outgroups
    for org in organisms:
        outgroups = set(organisms) - {org}
        job[org] = outgroups

    return job


def get_job(organisms: list[str], job_path: Optional[str | Path]) -> Dict[str, Set[str]]:
    """
    Resolve the LSE job definition.

    If a job file is provided, load it. Otherwise, use the default job where
    each organism uses all other organisms as outgroups.

    Args:
        organisms (list[str]): Known organism identifiers present in the dataset.
        job_path (str | Path | None): Path to the job description file, if any.

    Returns:
        Dict[str, Set[str]]: Mapping from query organism to its set of outgroups.
    """

    if job_path is None:
        return generate_default_job(organisms)

    job_path = Path(job_path)
    if not job_path.exists():
        raise FileNotFoundError(f"Job file not found: {job_path}")

    job = load_job_file(job_path)

    # Robustness: ignore organisms not present in the current dataset
    known = set(organisms)
    job = {q: {og for og in ogs if og in known} for q, ogs in job.items() if q in known}

    return job



# -----------------------------
# 3.  Processing
# -----------------------------

def add_org_columns(hits_df: pd.DataFrame, id_to_org: Dict[int, str]) -> pd.DataFrame:
    """
    Annotate BLAST hits with organism information.

    Adds the organism of the query and target proteins as new columns
    (org_q, org_t). This enables all downstream LSE filters to operate
    on organism-level relationships.
    """

    # Work on a copy to avoid modifying the original hits table
    df = hits_df.copy()

    # Map protein IDs to organism labels (NaN if ID is unknown)
    df["org_q"] = df["query"].map(id_to_org)
    df["org_t"] = df["subject"].map(id_to_org)

    return df



# -----------------------------
# 4.  Algorithm / Core
# -----------------------------


def early_stop_intra_hits_for_query_org(
    pre_grouped_hits,
    query_org: str,
    outgroups: Set[str],
    queries_by_org,
) -> pd.DataFrame:
    """
    Compute directed intra-genome hits for one query organism using early-stop.

    For each query protein from `query_org`, hits are scanned in best-first order
    (lowest evalue, then highest score). Intra-genome hits are kept until the
    first outgroup hit is encountered; after that, all remaining hits for that
    query are ignored.

    Args:
        pre_grouped_hits: Dict[int, list[tuple[int, str, float]]]
            Mapping query_id -> list of (subject_id, subject_org, score) already
            sorted best-first.
        query_org (str): Organism currently processed as query organism.
        outgroups (Set[str]): Organisms treated as outgroups for early-stop.
        queries_by_org: Dict[str, list[int]]
            Mapping organism -> list of query_ids belonging to that organism.

    Returns:
        pd.DataFrame: Directed intra-hits with columns ["query", "subject", "score"].
    """
    kept_rows = []

    for query_id in queries_by_org.get(query_org, []):
        hits = pre_grouped_hits.get(query_id, [])
            

        seen_subjects = set()

        for subject_id, org_t, score in hits:
            if subject_id in seen_subjects:
                continue
            seen_subjects.add(subject_id)

            # Outgroup-hit -> stop scanning this query immediately
            if org_t in outgroups:
                break

            # intra-gnome hit -> keep
            if org_t == query_org:
                kept_rows.append((query_id, subject_id, score))

    return pd.DataFrame(kept_rows, columns=["query", "subject", "score"])



def reciprocal_edges_from_directed_hits(directed_hits: pd.DataFrame) -> List[Tuple[int, int]]:
    """
    Build reciprocal (undirected) edges from directed intra-genome hits.

    A connection between two proteins is kept
    only if hits exist in both directions (p → q and q → p).
    
    Args:
        directed_hits (pd.DataFrame): Directed hits with columns ["query", "subject"].

    Returns:
        List[Tuple[int, int]]: Undirected edges represented as (min_id, max_id).
    """

    if directed_hits.empty:
        return []

    # Collect all directed (query, subject) pairs
    pairs = set(zip(directed_hits["query"].astype(int), directed_hits["subject"].astype(int)))

    edges: List[Tuple[int, int]] = []

    # Keep an edge only if both directions are present
    for a, b in pairs:
        if a == b:
            continue
        if (b, a) in pairs:
            u, v = (a, b) if a < b else (b, a)
            edges.append((u, v))

    # Deduplicate undirected edges
    edges = list(set(edges))

    return edges


def undirected_edge_allowed_by_presence(
        a: int,
        b: int, 
        presence: dict[int, list[int]]
        ) -> bool:
    """
    Check whether an undirected edge is allowed by strict presence.

    Used to filter undirected edges. An edge between
    two proteins is kept only if both directions are present in the presence
    map.
    """

    if a == b:
        return False
    
    # Require presence in both directions
    return (b in presence.get(a, [])) and (a in presence.get(b, []))


def connected_components(
    nodes: Iterable[int],
    edges: Iterable[Tuple[int, int]],
) -> List[List[int]]:
    """
    Compute connected components of an undirected graph.

    Used to turn protein nodes and undirected edges
    into groups (connected components). Each returned component is a sorted list
    of protein IDs.
    
    Args:
        nodes (Iterable[int]): All node IDs that should be part of the graph.
        edges (Iterable[Tuple[int, int]]): Undirected edges between node IDs.

    Returns:
        List[List[int]]: List of components; each component is a sorted list of IDs.
    """

    # Build adjacency list: node -> set(presence_map)
    adj: Dict[int, Set[int]] = {n: set() for n in nodes}
    for a, b in edges:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)

    visited: Set[int] = set()
    components: List[List[int]] = []

    # Traverse all nodes; start a graph search for each unvisited node
    for start in adj.keys():
        if start in visited:
            continue

        stack = [start]
        visited.add(start)
        comp: List[int] = []

        # Depth-first search to collect all nodes in this connected component
        while stack:
            v = stack.pop()
            comp.append(v)

            for nb in adj[v]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)

        # Sort component IDs for deterministic output
        components.append(sorted(comp))

    return components


def compute_lse_groups_for_org(
    id_to_org: Dict[int, str],
    query_org: str,
    outgroups: Set[str],
    presence_map: dict[int, list[int]],
    pre_grouped_hits,
    queries_by_org,
) -> List[List[int]]:
    """
    Compute LSE groups for a single query organism.

    Per-organism pipeline:
    (1) early-stop directed intra-hits,
    (2) reciprocal undirected edges,
    (3) strict presence filter,
    (4) connected components -> groups.

    Args:
        id_to_org (Dict[int, str]): Protein-ID -> organism mapping.
        query_org (str): Organism for which groups are computed.
        outgroups (Set[str]): Outgroup organisms for early-stop.
        presence_map (dict[int, list[int]]):
            Presence/adjacency constraints used by `undirected_edge_allowed_by_presence`.
        pre_grouped_hits: Dict[int, list[tuple[int, str, float]]]
            Mapping query_id -> list of (subject_id, subject_org, score), sorted best-first.
        queries_by_org: Dict[str, list[int]]
            Mapping organism -> list of query_ids.

    Returns:
        List[List[int]]: LSE groups as connected components (sorted ID lists).
    """

    # All protein nodes that belong to this organism
    nodes = sorted(pid for pid, org in id_to_org.items() if org == query_org)

    # Directed intra-hits after early-stop (per query protein)
    directed_intra = early_stop_intra_hits_for_query_org(
        pre_grouped_hits=pre_grouped_hits,
        query_org=query_org,
        outgroups=outgroups,
        queries_by_org=queries_by_org,
    )

    # Keep only reciprocal relationships (p->q and q->p) as undirected edges
    undirected_edges = reciprocal_edges_from_directed_hits(directed_intra)

    # Strict presence filter.
    undirected_edges = [
        (a, b) for (a, b) in undirected_edges
        if undirected_edge_allowed_by_presence(a, b, presence_map)
    ]


    # Connected components of the undirected graph are the final groups
    groups = connected_components(nodes=nodes, edges=undirected_edges)
    return groups




# -----------------------------
# 5.  Output
# -----------------------------

def order_groups_cpp_like(groups: List[List[int]]) -> List[List[int]]:
    """
    Order LSE groups according to the original COGlse output semantics.

    Groups with size >= 2 are listed first (in their discovery order),
    followed by singleton groups.
    
    Args:
        groups (List[List[int]]): LSE groups as lists of protein IDs.

    Returns:
        List[List[int]]: Reordered groups matching COGlse conventions.
    """

    multi = []
    single = []

    for g in groups:
        g_sorted = sorted(g)
        if len(g_sorted) >= 2:
            multi.append(g_sorted)
        else:
            single.append(g_sorted)

    return multi + single


def groups_to_df_lines(groups: List[List[int]]) -> pd.DataFrame:
    """
    Convert LSE groups into the final output DataFrame format.

    Groups are ordered according to COGlse conventions and serialized as
    comma-separated protein ID strings, one group per row.

    Args:
        groups (List[List[int]]): LSE groups as lists of protein IDs.

    Returns:
        pd.DataFrame: Output table with a single column "group".
    """

    groups = order_groups_cpp_like(groups)
    lines = [",".join(str(x) for x in g) for g in groups]

    return pd.DataFrame({"group": lines})



# -----------------------------
# 6.  Runner
# -----------------------------

def run_lse(
    hits_path: str | Path,
    prot2org_path: str | Path,
    presence_map: str | Path,
    output_path: str | Path | None = None,
    job_path: str | Path | None = None,
) -> pd.DataFrame:
    """
    lse subcommand entry point
    """

    hits = load_hits(hits_path)
    prot2org = load_prot2org(prot2org_path)
    id_to_org = build_id_to_org(prot2org)
    organisms = list_organisms(id_to_org)
    job = get_job(organisms, job_path)
    hits_with_org = add_org_columns(hits, id_to_org)
    presence_map = load_presence_map(presence_map)

    hits_with_org = hits_with_org.sort_values(
        ["query", "evalue", "score"],
        ascending=[True, True, False],
        kind="mergesort",
    )

    # pre_grouped_hits: query -> list[(subject, subject_org, score)]
    pre_grouped_hits = defaultdict(list)
    for row in hits_with_org.itertuples(index=False):
        pre_grouped_hits[row.query].append((row.subject, row.org_t, row.score))

    # queries_by_org: organism -> list[query_ids]
    queries_by_org = defaultdict(list)
    for q in pre_grouped_hits.keys():
        q_org = id_to_org.get(q)
        if q_org is not None:
            queries_by_org[q_org].append(q)

    all_groups: List[List[int]] = []
    for query_org, outgroups in job.items():
        groups = compute_lse_groups_for_org(
            id_to_org=id_to_org,
            query_org=query_org,
            outgroups=outgroups,
            presence_map=presence_map,
            pre_grouped_hits=pre_grouped_hits,
            queries_by_org=queries_by_org,
        )
        all_groups.extend(groups)


    lse_df = groups_to_df_lines(all_groups)

    # optional: write file (no Header, no index)
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        lse_df["group"].to_csv(output_path, index=False, header=False)

    return lse_df

def run_lse_in_memory(
    hits,
    prot2org ,
    presence_map: dict[int, list[int]],
    output_path: str | Path | None = None,
    job_path: str | Path | None = None,
) -> pd.DataFrame:
    
    """

    Run LSE from in-memory tables.

    This is the pipeline entry point used by the Python run pipeline.
    """

    id_to_org = build_id_to_org(prot2org)
    organisms = list_organisms(id_to_org)
    job = get_job(organisms, job_path)
    hits_with_org = add_org_columns(hits, id_to_org)

    # global stable sort
    hits_with_org = hits_with_org.sort_values(
        ["query", "evalue", "score"],
        ascending = [True, True, False],
        kind = "mergesort",
    )

    # Pregroup hits per query
    pre_grouped_hits = defaultdict(list)
    for row in hits_with_org.itertuples(index=False):
        pre_grouped_hits[row.query].append((row.subject, row.org_t, row.score))


    queries_by_org = defaultdict(list)
    for query_id in pre_grouped_hits.keys():
        org = id_to_org.get(query_id)
        if org is not None:
            queries_by_org[org].append(query_id)


    logger.info(
        "Setup complete: hits=%d, proteins=%d, organisms=%d, job_queries=%d",
        len(hits),
        len(id_to_org),
        len(organisms),
        len(job),
    )

    all_groups: List[List[int]] = []

    for query_org, outgroups in job.items():
        logger.debug("%s: start (outgroups=%d)", query_org, len(outgroups))

        groups = compute_lse_groups_for_org(
            id_to_org=id_to_org,
            query_org=query_org,
            outgroups=outgroups,
            presence_map=presence_map,  
            pre_grouped_hits=pre_grouped_hits,
            queries_by_org=queries_by_org,
        )
        
        n_groups = len(groups)
        n_members = sum(len(g) for g in groups)
        logger.debug("%s: done (clusters=%d, proteins=%d)", query_org, n_groups, n_members)

        all_groups.extend(groups)


    lse_df = groups_to_df_lines(all_groups)
    logger.info("Output built: total_groups=%d", len(lse_df))

    # optional: write file (no Header, no index)
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        lse_df["group"].to_csv(output_path, index=False, header=False)
        logger.info("Output written: %s", output_path)

    return lse_df
    