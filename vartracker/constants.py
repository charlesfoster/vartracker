"""
Constants and mappings for SARS-CoV-2 analysis.

Contains genome coordinates, gene mappings, and other SARS-CoV-2 specific data.
"""

import re
from bisect import bisect_right
from typing import List, cast

# NSP (non-structural protein) coordinates and mappings for SARS-CoV-2
NSPS = {
    "genes": [
        "nsp1",
        "nsp2",
        "nsp3",
        "nsp4",
        "nsp5",
        "nsp6",
        "nsp7",
        "nsp8",
        "nsp9",
        "nsp10",
        "nsp11",
        "nsp12",
        "nsp13",
        "nsp14",
        "nsp15",
        "nsp16",
    ],
    "aa_start": [
        1,
        181,
        819,
        2764,
        3264,
        3570,
        3860,
        3943,
        4141,
        4254,
        4393,
        4393,
        5325,
        5926,
        6453,
        6799,
    ],
    "nt_start": [
        266,
        806,
        2720,
        8555,
        10055,
        10973,
        11843,
        12092,
        12686,
        13025,
        13442,
        13442,
        16237,
        18040,
        19621,
        20659,
    ],
    "nt_end": [
        805,
        2719,
        8554,
        10054,
        10972,
        11842,
        12091,
        12685,
        13024,
        13441,
        13480,
        16236,
        18039,
        19620,
        20658,
        21552,
    ],
    "product": [
        "nsp1",
        "nsp2",
        "nsp3_PLpro",
        "nsp4_TM",
        "nsp5_3CLpro",
        "nsp6_TM",
        "nsp7",
        "nsp8",
        "nsp9",
        "nsp10_CysHis",
        "nsp11",
        "nsp12_RdRp",
        "nsp13_ZBD",
        "nsp14_exonuclease",
        "nsp15_NendoU",
        "nsp16_0MT",
    ],
}

# Reference gene lengths for SARS-CoV-2
REF_GENE_LENGTHS = {
    "5' UTR": 265,
    "ORF1ab": 21291,
    "S": 3822,
    "ORF3a": 828,
    "ORF3c": 126,
    "E": 228,
    "M": 669,
    "ORF6": 186,
    "ORF7a": 366,
    "ORF7b": 132,
    "ORF8": 366,
    "N": 1260,
    "ORF9b": 294,
    "ORF10": 117,
    "3' UTR": 228,
    "INTERGENIC": 1,
}

# Calculate NSP lengths
NSP_LENGTHS = {}
nsp_products = cast(List[str], NSPS["product"])
nsp_nt_start = cast(List[int], NSPS["nt_start"])
nsp_nt_end = cast(List[int], NSPS["nt_end"])

for nsp_name, start, end in zip(nsp_products, nsp_nt_start, nsp_nt_end):
    NSP_LENGTHS[nsp_name] = end - start + 1


def bcf_orf1ab_to_nsp(mutation):
    """
    Convert ORF1ab mutation to NSP-specific notation.

    Args:
        mutation (str): Mutation string from bcftools csq

    Returns:
        str: NSP-specific mutation notation

    Raises:
        ValueError: If mutation format is invalid
    """
    # Parse out important bits from mutation
    no_gene = re.sub(".*:", "", mutation)
    parsed = re.findall(r"([A-Z]+)", no_gene)

    if len(parsed) == 2:
        ref = parsed[0]
        alt = parsed[1]
    elif len(parsed) == 1:
        ref = ""
        alt = parsed[0]
    else:
        raise ValueError(f"Invalid mutation format: {mutation}")

    pos_match = re.search(r"(\d+)", no_gene)
    if not pos_match:
        raise ValueError(f"No position found in mutation: {mutation}")

    pos = int(pos_match.group())

    # Find corresponding nsp
    idx = bisect_right(NSPS["aa_start"], pos) - 1
    if idx < 0 or idx >= len(NSPS["product"]):
        raise ValueError(f"Position {pos} out of range for NSP mapping")

    nsp = NSPS["product"][idx]
    new_pos = pos - NSPS["aa_start"][idx] + 1
    new_variant = f"{nsp}:{ref}{new_pos}{alt}"

    return new_variant


def scrape_orf1ab_to_nsp(nt_mutation, aa_mutation):
    """
    Convert ORF1ab mutations to NSP-specific notation for scraped data.

    Args:
        nt_mutation (str): Nucleotide mutation
        aa_mutation (str): Amino acid mutation

    Returns:
        tuple: (nsp, nt_variant, aa_variant)
    """
    # Process AA mutation
    parsed = re.findall(r"([A-Z]+)", aa_mutation)
    if len(parsed) == 2:
        ref = parsed[0]
        alt = parsed[1]
    elif len(parsed) == 1:
        ref = ""
        alt = parsed[0]
    else:
        ref = ""
        alt = ""

    pos_match = re.search(r"(\d+)", aa_mutation)
    pos = int(pos_match.group()) if pos_match else 0

    # Find corresponding nsp for AA
    if 4393 <= pos <= 4405:
        nsp = "nsp11/nsp12_RdRp"
        new_pos = pos - NSPS["aa_start"][10] + 1
    else:
        idx = bisect_right(NSPS["aa_start"], pos) - 1
        nsp = NSPS["genes"][idx] if 0 <= idx < len(NSPS["genes"]) else "unknown"
        new_pos = (
            pos - NSPS["aa_start"][idx] + 1 if 0 <= idx < len(NSPS["aa_start"]) else pos
        )

    new_aa_variant = f"{ref}{new_pos}{alt}"

    # Process NT mutation
    parsed_nt = re.findall(r"([A-Z]+)", nt_mutation)
    if len(parsed_nt) == 2:
        ref_nt = parsed_nt[0]
        alt_nt = parsed_nt[1]
    elif len(parsed_nt) == 1:
        ref_nt = ""
        alt_nt = parsed_nt[0]
    else:
        ref_nt = ""
        alt_nt = ""

    pos_nt_match = re.search(r"(\d+)", nt_mutation)
    pos_nt = int(pos_nt_match.group()) if pos_nt_match else 0

    # Find corresponding nsp for NT
    if 13442 <= pos_nt <= 13480:
        new_pos_nt = pos_nt - NSPS["nt_start"][10] + 1
    else:
        idx_nt = bisect_right(NSPS["nt_start"], pos_nt) - 1
        new_pos_nt = (
            pos_nt - NSPS["nt_start"][idx_nt] + 1
            if 0 <= idx_nt < len(NSPS["nt_start"])
            else pos_nt
        )

    new_nt_variant = f"{ref_nt}{new_pos_nt}{alt_nt}"

    return (nsp, new_nt_variant, new_aa_variant)


def reformat_csq_notation(gene, string):
    """
    Reformat bcftools csq notation to standard format.

    Args:
        gene (str): Gene name
        string (str): CSQ notation string

    Returns:
        tuple: (reformatted_string, nsp_info)
    """
    if not string:
        return ("", "")

    splitter = string.rfind(">")
    if splitter == -1:
        return (string, "")

    ref = re.sub("[0-9]", "", string[:splitter])
    pos = re.sub("[A-Za-z]", "", string[:splitter])
    alt = re.sub("[0-9]", "", string[splitter + 1 :])
    reformatted = ref + pos + alt

    if gene == "ORF1ab":
        try:
            nsp_info = bcf_orf1ab_to_nsp(reformatted)
        except Exception:
            nsp_info = "ERROR"
    else:
        nsp_info = ""

    return (reformatted, nsp_info)
