import re

def extract_organism(header: str) -> str:
    """Extract the organism name from a Fasta definition line (line starting with '>')

    Args:
        header (str): full definition line

    Raises:
        ValueError: If multiple matches are found for the organism modifier in the header '[organism=...]'.

    Returns:
        str: The extracted organism name. If no explicit organism modifier is found, it attempts to extract a bare genus/species from the rightmost brackets without a modifier key. If no organism information can be extracted, it returns "Unknown".
    """

    modifier_keys = {
        "organism", "strain", "chromosome", "sex", "tissue-type", "moltype", "altitude",
        "bio-material", "breed", "cell-line", "cell-type", "clone", "collected-by",
        "collection-date", "cultivar", "culture-collection", "dev-stage", "ecotype",
        "endogenous-virus-name", "fwd-pcr-primer-name", "fwd-pcr-primer-seq", "genotype",
        "geo-loc-name", "haplogroup", "haplotype", "host", "isolate", "isolation-source",
        "lab-host", "lat-lon", "linkage-group", "map", "mating-type", "note", "plasmid-name",
        "plastid-name", "rev-pcr-primer-name", "rev-pcr-primer-seq", "segment", "serotype",
        "serovar", "specimen-voucher", "sub-species", "variety", "gcode", "mgcode",
        "environmental-sample", "germline", "metagenomic", "rearranged", "transgenic"
    }

    # Check if any explicit [organism=...] is present
    explicit = re.findall(r'\[organism=([^\]]+)\]', header, re.IGNORECASE)
    if explicit:
        if len(explicit) > 1:
            raise ValueError(f"Multiple matches found for organism in header: {header}")
        elif len(explicit) == 1:
            return explicit[0].strip()

    # If no explicit [organism=...] is found, we look at all other brackets
    all_brackets = list(re.finditer(r'\[([^\]]+)\]', header))

    # block_start tracks the leftmost position of any confirmed modifier block
    block_start = len(header)

    # Iterate from the rightmost bracket to the left, since NCBI modifiers are supposed to be at the end of the header
    for match in reversed(all_brackets):
        content = match.group(1).strip()
        key = content.split("=")[0].strip().lower()

        is_keyed = "=" in content and key in modifier_keys
        is_bare = "=" not in content

        # Check that everything after this bracket is whitespace or in the block
        text_after = header[match.end():block_start].strip()
        is_at_tail = text_after == "" or all_brackets_in_range(header, match.end(), block_start)

        if is_keyed and is_at_tail:
            # This bracket is part of a modifier block with a different key than organism, so we exclude it from the "is_at_tail" check for the next iteration.
            block_start = match.start()

        elif is_bare and is_at_tail:
            # This bracket is a bare genus/species at the end of the header, so we extract it and return. We take at most the first two words to avoid strains or other modifiers that might be included in the same bracket.
            genus_species = content.split(" ", 2)
            return genus_species[0] if len(genus_species) == 1 else " ".join(genus_species[:2])
        else:
            # If we encounter a bracket that doesn't fit the modifier block, we stop looking further left.
            break  

    return "Unknown"


def all_brackets_in_range(header: str, start: int, end: int) -> bool:
    """Check that the substring header[start:end] contains only whitespace and bracket spans that are already known.

    Args:
        header (str): The full string to check within.
        start (int): The starting index of the substring to check.
        end (int): The ending index of the substring to check.

    Returns:
        bool: True if the substring contains only whitespace and known bracket spans, False otherwise.
    """
    segment = header[start:end]
    cleaned = re.sub(r'\[[^\]]+\]', '', segment)
    return cleaned.strip() == ""


headers = [
    ">NP_219502.1 hypothetical protein CT_875 [Chlamydia trachomatis D/UW-3/CX]",
    ">NP_415603.2 DUF2655 domain-containing protein YceQ [Escherichia coli str. K-12 substr. MG1655]",
    ">NP_415609.1 3-oxoacyl-[acyl carrier protein] synthase 3 [Escherichia coli str. K-12 substr. MG1655]",
    ">WP_000009173.1 MULTISPECIES: YeiH family protein [Streptococcus]",
    ">WP_000011013.1 DEAD/DEAH box helicase [Streptococcus pneumoniae]",
    ">YP_008718.1 hypothetical protein [Streptococcus pneumoniae] [strain=R6]",
    ">YP_008718.1 hypothetical protein [organism=Streptococcus pneumoniae] [strain=R6]",
    ">YP_008718.1 hypothetical protein [organism=Streptococcus pneumoniae] [strain=R6] [chromosome=1]",
    ">YP_008718.1 hypothetical protein [organism=Streptococcus] [strain=R6]",
    ">NP_219502.1 hypothetical protein CT_875",
    ">NP_415609.1 3-oxoacyl-[acyl carrier protein] synthase 3",
    ">YP_008718.1 hypothetical protein [Streptococcus pneumoniae R6] [chromosome=1]",
    ">NP_414716.1 ditrans,polycis-undecaprenyl-diphosphate synthase [(2E,6E)-farnesyl-diphosphate specific] [Escherichia coli str. K-12 substr. MG1655]"#,
    #">YP_008718.1 hypothetical protein [organism=Streptococcus] [organism=R6]"
]


if __name__ == "__main__":
    
    for header in headers:
        organism = extract_organism(header)
        print(f"Header: {header}\nExtracted organism: {organism}\n\n")
