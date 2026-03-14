import pytest
from COGprocessinput.extract_org import extract_organism

valid_test_headers = [
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
    ">NP_414716.1 ditrans,polycis-undecaprenyl-diphosphate synthase [(2E,6E)-farnesyl-diphosphate specific] [Escherichia coli str. K-12 substr. MG1655]"
]

invalid_test_headers = [
    ">NP_219502.1 hypothetical protein CT_875 [organism=Chlamydia trachomatis D/UW-3/CX] [organism=Escherichia coli]"
]

expected_organisms = [
    "Chlamydia trachomatis",
    "Escherichia coli",
    "Escherichia coli",
    "Streptococcus",
    "Streptococcus pneumoniae",
    "Streptococcus pneumoniae",
    "Streptococcus pneumoniae",
    "Streptococcus pneumoniae",
    "Streptococcus",
    "Unknown",
    "Unknown",
    "Streptococcus pneumoniae",
    "Escherichia coli"
]

def test_extract_organism_valid():
    for header, expected in zip(valid_test_headers, expected_organisms):
        assert extract_organism(header) == expected, f"Header: {header}\nExpected: {expected}\nGot: {extract_organism(header)}\n\n"

def test_extract_organism_invalid():
    with pytest.raises(ValueError):
        for header in invalid_test_headers:
            extract_organism(header)