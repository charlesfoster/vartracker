#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: github.com/charlesfoster
"""

import os
import shutil
import re
from random import randrange
import subprocess
import argparse
import sys
import pandas as pd
from cyvcf2 import VCF, Writer
import matplotlib.pyplot as plt

plt.rcdefaults()
import string
import seaborn as sns
from bisect import bisect_right
import numpy as np
from itertools import dropwhile
from argparse_formatter import FlexiFormatter

# logo
logo = '''
██    ██  █████  ██████  ████████ ██████   █████   ██████ ██   ██ ███████ ██████  
██    ██ ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██ 
██    ██ ███████ ██████     ██    ██████  ███████ ██      █████   █████   ██████  
 ██  ██  ██   ██ ██   ██    ██    ██   ██ ██   ██ ██      ██  ██  ██      ██   ██ 
  ████   ██   ██ ██   ██    ██    ██   ██ ██   ██  ██████ ██   ██ ███████ ██   ██ 
                                                                                  
                                                                                  
'''
# key variables
thisdir = os.path.abspath(os.path.dirname(__file__))

nsps = {
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
        "nsp3",
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


# functions
def rindex(lst, item):
    def index_ne(x):
        return lst[x] != item

    try:
        return next(dropwhile(index_ne, reversed(range(len(lst)))))
    except StopIteration:
        raise ValueError("rindex(lst, item): item not in list")


def format_vcf(
    vcf,
    sample,
    tempdir,
    rename_chrom,
    chrom_map,
    min_snv_freq,
    min_indel_freq,
    reference,
    annotation,
    debug,
):
    out = os.path.join(tempdir, os.path.basename(vcf))
    out = re.sub(".vcf.gz", ".cyvcf2.vcf.gz", out)
    stem = sample
    log = os.path.join(tempdir, stem + ".log")
    csq_file = re.sub(".cyvcf2.vcf.gz", ".csq.vcf.gz", out)
    vcf_mod = VCF(vcf, strict_gt=True)
    vcf_mod.add_format_to_header(
        {"ID": "DP", "Description": "Raw Depth", "Type": "Integer", "Number": "1"}
    )
    vcf_mod.add_format_to_header(
        {"ID": "AF", "Description": "Allele Frequency", "Type": "Float", "Number": "1"}
    )
    w = Writer(out, vcf_mod, mode="wz")
    variants = {}
    for v in vcf_mod:
        if str(v.POS) in variants.keys():
            if variants[str(v.POS)].INFO["AF"] < v.INFO["AF"]:
                del variants[str(v.POS)]
            else:
                continue
        # curate variant
        v.ID = v.REF + str(v.POS) + v.ALT[0]
        v.set_format("DP", np.array([v.INFO["DP"]]))
        v.set_format("AF", np.array([v.INFO["AF"]]))
        v.genotypes = [[1, True]]
        # add to variants dict
        variants[str(v.POS)] = v
    for v in variants.items():
        w.write_record(v[1])
    w.close()
    # bgzip and index cyout
    cmd = f"bcftools index -f {out}"
    os.system(cmd)
    # run csq before merging so that AF filtering is more appropriate
    cmd = 'bcftools view -i \'(INFO/AF >= {0} & INFO/TYPE == "SNP") | (INFO/AF >= {1} & INFO/TYPE == "INDEL")\' {2} | bcftools csq -f {3} -g {4} --force -Oz -o {5}; tabix -f -p vcf {5}'.format(
        min_snv_freq, min_indel_freq, out, reference, annotation, csq_file
    )
    with open(log, "a", encoding="utf-8") as err:
        subprocess.Popen(
            cmd,
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=err,
            universal_newlines=True,
        ).communicate()
    if debug:
        print(repr(cmd))


def merge_consequences(tempdir, csq_file, sample_names, debug):
    cmd = "bcftools merge -0 -Ov -l {0} | bcftools reheader -s {1} > {2}".format(
        os.path.join(tempdir, "vcf_list.txt"), sample_names, csq_file
    )
    if debug:
        print(repr(cmd))
    subprocess.Popen(
        cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    ).communicate()


def reformat_csq_notation(gene, string):
    splitter = string.rfind(">")
    ref = re.sub("[0-9]", "", string[:splitter])
    pos = re.sub("[A-Za-z]", "", string[:splitter])
    alt = re.sub("[0-9]", "", string[splitter + 1 :])
    reformatted = ref + pos + alt
    if gene == "ORF1ab":
        try:
            nsp_info = bcf_orf1ab_to_nsp(reformatted)
        except:
            nsp_info = "ERROR"
    else:
        nsp_info = ""
    return (reformatted, nsp_info)


def bcf_orf1ab_to_nsp(mutation):
    # parse out important bits from mutation
    no_gene = re.sub(".*:", "", mutation)
    parsed = re.findall(r"([A-Z]+)", no_gene)
    if len(parsed) == 2:
        ref = parsed[0]
        alt = parsed[1]
    else:
        ref = ""
        alt = parsed[0]
    pos = re.search(r"(\d)+", no_gene).group()
    # find corresponding nsp
    idx = bisect_right(nsps["aa_start"], int(pos)) - 1
    nsp = nsps["product"][idx]
    new_pos = int(pos) - nsps["aa_start"][idx] + 1
    new_variant = f"{nsp}:{ref}{new_pos}{alt}"
    return new_variant


def scrape_orf1ab_to_nsp(nt_mutation, aa_mutation):
    ## AA
    # parse out important bits from mutation
    parsed = re.findall(r"([A-Z]+)", aa_mutation)
    if len(parsed) == 2:
        ref = parsed[0]
        alt = parsed[1]
    else:
        ref = ""
        alt = parsed[0]
    pos = re.search(r"(\d)+", aa_mutation).group()
    # find corresponding nsp
    if int(pos) >= 4393 and int(pos) <= 4405:
        nsp = "nsp11/nsp12_RdRp"
        new_pos = int(pos) - nsps["aa_start"][10] + 1
    else:
        idx = bisect_right(nsps["nt_start"], int(pos)) - 1
        nsp = nsps["genes"][idx]
        new_pos = int(pos) - nsps["nt_start"][idx] + 1
    idx = bisect_right(nsps["aa_start"], int(pos)) - 1
    nsp = nsps["genes"][idx]
    new_pos = int(pos) - nsps["aa_start"][idx] + 1
    new_aa_variant = f"{ref}{new_pos}{alt}"
    ## NT
    # parse out important bits from mutation
    parsed = re.findall(r"([A-Z]+)", nt_mutation)
    if len(parsed) == 2:
        ref = parsed[0]
        alt = parsed[1]
    else:
        ref = ""
        alt = parsed[0]
    pos = re.search(r"(\d)+", nt_mutation).group()
    # find corresponding nsp
    if int(pos) >= 13442 and int(pos) <= 13480:
        nsp = "nsp11/nsp12_RdRp"
        new_pos = int(pos) - nsps["nt_start"][10] + 1
    else:
        idx = bisect_right(nsps["nt_start"], int(pos)) - 1
        nsp = nsps["genes"][idx]
        new_pos = int(pos) - nsps["nt_start"][idx] + 1
    new_nt_variant = f"{ref}{new_pos}{alt}"
    return (nsp, new_nt_variant, new_aa_variant)


def calculate_variant_site_depths(cov_df, v, samples):
    start = v.start + 2  # account for 0-based coords and wanting the first alt base
    if v.is_indel and not v.is_deletion:
        end = v.end + 2
    else:
        end = v.end
    if v.var_type == "snp":
        site_depths = list(cov_df.loc[cov_df["pos"] == start]["depth"])
    else:
        site_depths = (
            cov_df[cov_df["pos"].between(start, end, inclusive="both")]
            .groupby("sample")["depth"]
            .mean()
            .reset_index()
        )
        site_depths["order"] = pd.Categorical(
            site_depths["sample"], categories=samples, ordered=True
        )
        site_depths = [
            int(round(x, 0)) for x in list(site_depths.sort_values("order")["depth"])
        ]
    window_depth = (
        cov_df[cov_df["pos"].between(start - 10, end + 10, inclusive="neither")]
        .groupby("sample")["depth"]
        .mean()
        .reset_index()
    )
    window_depth["order"] = pd.Categorical(
        window_depth["sample"], categories=samples, ordered=True
    )
    window_depth = [
        int(round(x, 0)) for x in list(window_depth.sort_values("order")["depth"])
    ]
    variant_depths = [
        (x[0] * 0) - 1 if x[0] < 0 else x[0] for x in v.format("DP").tolist()
    ]
    variant_qc = []
    for x, y, z in zip(variant_depths, site_depths, window_depth):
        if x >= 0:
            variant_qc.append("P")
        elif x < 0 and y >= 10:
            variant_qc.append("P")
        elif v.is_indel and x < 0 and y < 10 and z >= 10:
            variant_qc.append("P")
        elif v.var_type == "snp" and x < 0 and y < 10:
            variant_qc.append("F")
        elif x < 0 and y < 10 and z < 10:
            variant_qc.append("F")
    variant_depths = ["M" if x == -1 else x for x in variant_depths]
    result = {
        "site_depth": [str(x) for x in site_depths],
        "window_depth": [str(x) for x in window_depth],
        "variant_depth": [str(x) for x in variant_depths],
        "variant_qc": [str(x) for x in variant_qc],
    }
    return result


class AminoAcidChange:
    def __init__(self, aa_change):
        amino_acids = {
            "I": ["Hydrophobic", "Aliphatic", "Uncharged"],
            "L": ["Hydrophobic", "Aliphatic", "Uncharged"],
            "V": ["Hydrophobic", "Small", "Aliphatic", "Uncharged"],
            "C": ["Hydrophobic", "Small", "Uncharged"],
            "A": ["Hydrophobic", "Small", "Tiny", "Uncharged"],
            "G": ["Hydrophobic", "Small", "Tiny", "Uncharged"],
            "M": ["Hydrophobic", "Uncharged"],
            "F": ["Hydrophobic", "Aromatic", "Uncharged"],
            "Y": ["Hydrophobic", "Polar", "Aromatic", "Uncharged"],
            "W": ["Hydrophobic", "Polar", "Aromatic", "Uncharged"],
            "H": ["Hydrophobic", "Polar", "Aromatic", "Positive", "Charged"],
            "K": ["Hydrophobic", "Polar", "Positive", "Charged"],
            "R": ["Polar", "Positive", "Charged"],
            "E": ["Polar", "Negative", "Charged"],
            "Q": ["Polar", "Uncharged"],
            "D": ["Polar", "Proline", "Negative", "Charged"],
            "N": ["Polar", "Proline", "Uncharged"],
            "S": ["Polar", "Proline", "Tiny", "Uncharged"],
            "T": ["Hydrophobic", "Polar", "Proline", "Uncharged"],
            "P": ["Proline", "Uncharged"],
            "B": ["Polar", "Uncharged"],
            "Z": ["Polar", "Uncharged"],
            "X": [
                "Hydrophobic",
                "Polar",
                "Tiny",
                "Aliphatic",
                "Aromatic",
                "Positive",
                "Negative",
                "Charged",
                "Uncharged",
            ],
            "-": [
                "Hydrophobic",
                "Polar",
                "Tiny",
                "Aliphatic",
                "Aromatic",
                "Positive",
                "Negative",
                "Charged",
                "Uncharged",
            ],
        }

        properties = {
            "Hydrophobic": [
                "I",
                "L",
                "V",
                "C",
                "A",
                "G",
                "M",
                "F",
                "Y",
                "W",
                "H",
                "K",
                "T",
                "X",
                "-",
            ],
            "Polar": [
                "Y",
                "W",
                "H",
                "K",
                "R",
                "E",
                "Q",
                "D",
                "N",
                "S",
                "T",
                "B",
                "Z",
                "X",
                "-",
            ],
            "Small": [
                "V",
                "C",
                "A",
                "G",
            ],
            "Proline": ["D", "N", "S", "T", "P"],
            "Tiny": ["A", "G", "S", "X", "-"],
            "Aliphatic": ["I", "L", "V", "X", "-"],
            "Aromatic": ["F", "Y", "W", "H", "X", "-"],
            "Positive": ["H", "K", "R", "X", "-"],
            "Negative": ["E", "D", "X", "-"],
            "Charged": ["H", "K", "R", "X", "E", "D", "-"],
            "Uncharged": [
                "I",
                "L",
                "V",
                "C",
                "A",
                "G",
                "M",
                "F",
                "Y",
                "W",
                "Q",
                "N",
                "S",
                "T",
                "P",
                "B",
                "Z",
                "X",
                "-",
            ],
        }

        weights = {
            "A": 71.037114,
            "R": 156.101111,
            "N": 114.042927,
            "D": 115.026943,
            "C": 103.009185,
            "E": 129.042593,
            "Q": 128.058578,
            "G": 57.021464,
            "H": 137.058912,
            "I": 113.084064,
            "L": 113.084064,
            "K": 128.094963,
            "M": 131.040485,
            "F": 147.068414,
            "P": 97.052764,
            "S": 87.032028,
            "T": 101.047679,
            "U": 150.95363,
            "W": 186.079313,
            "Y": 163.06332,
            "V": 99.068414,
            "X": 0,
            "B": 113.084064,
            "Z": 0,
        }

        try:
            aa1 = re.match(r"([a-z]+)([0-9]+)([a-z+])", aa_change, re.I).group(1)
            aa2 = re.match(r"([a-z]+)([0-9]+)([a-z+])", aa_change, re.I).group(3)
            self.change = aa_change
            self.aa1 = aa1
            self.aa2 = aa2
            self.aa1_total_properties = set(
                [
                    item
                    for sublist in [amino_acids.get(x) for x in aa1]
                    for item in sublist
                ]
            )
            self.aa2_total_properties = set(
                [
                    item
                    for sublist in [amino_acids.get(x) for x in aa2]
                    for item in sublist
                ]
            )
            self.aa1_unique_properties = list(
                set(self.aa1_total_properties - self.aa2_total_properties)
            )
            self.aa2_unique_properties = list(
                set(self.aa2_total_properties - self.aa1_total_properties)
            )
            self.aa1_weight = float("{:.3f}".format(sum([weights.get(x) for x in aa1])))
            self.aa2_weight = float("{:.3f}".format(sum([weights.get(x) for x in aa2])))
            self.weight_difference = round(self.aa1_weight - self.aa2_weight, 3) * -1
        except:
            self.change = aa_change
            self.aa1 = ""
            self.aa2 = ""
            self.aa1_total_properties = ""
            self.aa2_total_properties = ""
            self.aa1_unique_properties = ""
            self.aa2_unique_properties = ""
            self.aa1_weight = ""
            self.aa2_weight = ""
            self.weight_difference = float(0)

    def __repr__(self):
        print(f"Information for the following change: {self.change}")
        return str(self.__dict__)


def process_vcf(VCF_IN, covs):
    results = []
    vcf = VCF(VCF_IN)
    samples = vcf.samples
    cov_list = []
    total_cov_list = []
    for i in range(len(samples)):
        cov = pd.read_csv(covs[i], sep="\t", names=["ref", "pos", "depth"])
        total_cov = "{:.2f}".format(((sum(cov.depth > 10) / cov.shape[0]) * 100))
        total_cov_list.append(total_cov)
        cov["sample"] = samples[i]
        cov_list.append(cov)
    cov_df = pd.concat(cov_list).set_index("sample")
    for v in vcf:
        info = dict(list(v.INFO))
        allele_freqs = [
            "{:.3f}".format(x[0]).replace("nan", ".") for x in v.format("AF").tolist()
        ]
        presence_absence = ["N" if x == "." else "Y" for x in allele_freqs]
        variant_status = "new" if allele_freqs[0] == "." else "original"
        # calculate persistence of new variants
        if allele_freqs[0] != "." and allele_freqs[-1] == ".":
            persistent_status = "original_lost"
        elif allele_freqs[0] != "." and allele_freqs[-1] != ".":
            persistent_status = "original_retained"
        elif allele_freqs[0] == "." and allele_freqs[-1] != ".":
            persistent_status = "new_persistent"
        elif allele_freqs[0] == "." and allele_freqs[-1] == ".":
            persistent_status = "new_transient"
        depths_qc = calculate_variant_site_depths(cov_df, v, samples)
        overall_variant_qc = "FAIL" if "F" in depths_qc["variant_qc"] else "PASS"
        first_appearance = samples[presence_absence.index("Y")]
        last_appearance = samples[rindex(presence_absence, "Y")]
        if "BCSQ" in info.keys():
            anno = v.INFO["BCSQ"].split(",")[0].split("|")
            if len(anno) == 1 and anno[0].startswith("@") == True:
                result = {
                    "chrom": v.CHROM,
                    "start": v.start + 1,
                    "end": v.end,
                    "gene": anno[0],
                    "ref": v.REF,
                    "alt": v.ALT[0],
                    "variant": v.REF + str(v.POS) + v.ALT[0],
                    "amino_acid_consequence": anno[0],
                    "nsp_aa_change": anno[0],
                    "bcsq_nt_notation": anno[0],
                    "bcsq_aa_notation": anno[0],
                    "type_of_variant": v.var_type,
                    "type_of_change": anno[0],
                    "variant_status": variant_status,
                    "persistence_status": persistent_status,
                    "presence_absence": " / ".join(presence_absence),
                    "first_appearance": first_appearance,
                    "last_appearance": last_appearance,
                    "overall_variant_qc": overall_variant_qc,
                    "per_sample_variant_qc": " / ".join(depths_qc["variant_qc"]),
                    "aa1_total_properties": anno[0],
                    "aa2_total_properties": anno[0],
                    "aa1_unique_properties": anno[0],
                    "aa2_unique_properties": anno[0],
                    "aa1_weight": anno[0],
                    "aa2_weight": anno[0],
                    "weight_difference": anno[0],
                    "alt_freq": " / ".join(allele_freqs),
                    "variant_depth": " / ".join(depths_qc["variant_depth"]),
                    "variant_site_depth": " / ".join(depths_qc["site_depth"]),
                    "variant_window_depth": " / ".join(depths_qc["window_depth"]),
                    "samples": " / ".join(samples),
                    "total_genome_coverage": " / ".join(total_cov_list),
                }
            else:
                reformatted_nt = (
                    reformat_csq_notation(None, anno[6])
                    if len(anno) > 5
                    else [v.REF + str(v.POS) + v.ALT[0], ""]
                )
                reformatted_aa = (
                    reformat_csq_notation(anno[1], anno[5])
                    if len(anno) > 5
                    else [anno[1] + ":" + anno[0], ""]
                )
                aa_exploration = AminoAcidChange(reformatted_aa[0].replace("*", ""))

                result = {
                    "chrom": v.CHROM,
                    "start": v.start + 1,
                    "end": v.end,
                    "gene": anno[1],
                    "ref": v.REF,
                    "alt": v.ALT[0],
                    "variant": v.REF + str(v.POS) + v.ALT[0],
                    "amino_acid_consequence": reformatted_aa[0],
                    "nsp_aa_change": reformatted_aa[1],
                    "bcsq_nt_notation": anno[6] if len(anno) > 5 else "",
                    "bcsq_aa_notation": anno[5] if len(anno) > 5 else "",
                    "type_of_variant": v.var_type,
                    "type_of_change": anno[0],
                    "variant_status": variant_status,
                    "persistence_status": persistent_status,
                    "presence_absence": " / ".join(presence_absence),
                    "first_appearance": first_appearance,
                    "last_appearance": last_appearance,
                    "overall_variant_qc": overall_variant_qc,
                    "per_sample_variant_qc": " / ".join(depths_qc["variant_qc"]),
                    "aa1_total_properties": ";".join(
                        aa_exploration.aa1_total_properties
                    ),
                    "aa2_total_properties": ";".join(
                        aa_exploration.aa2_total_properties
                    ),
                    "aa1_unique_properties": ";".join(
                        aa_exploration.aa1_unique_properties
                    ),
                    "aa2_unique_properties": ";".join(
                        aa_exploration.aa2_unique_properties
                    ),
                    "aa1_weight": aa_exploration.aa1_weight,
                    "aa2_weight": aa_exploration.aa2_weight,
                    "weight_difference": aa_exploration.weight_difference,
                    "alt_freq": " / ".join(allele_freqs),
                    "variant_depth": " / ".join(depths_qc["variant_depth"]),
                    "variant_site_depth": " / ".join(depths_qc["site_depth"]),
                    "variant_window_depth": " / ".join(depths_qc["window_depth"]),
                    "samples": " / ".join(samples),
                    "total_genome_coverage": " / ".join(total_cov_list),
                }
            results.append(result)
        else:
            # continue
            if v.start + 1 < 266:
                gene = "5' UTR"
            elif v.start + 1 > 29674:
                gene = "3' UTR"
            else:
                gene = "INTERGENIC"
            result = {
                "chrom": v.CHROM,
                "start": v.start + 1,
                "end": v.end,
                "gene": gene,
                "ref": v.REF,
                "alt": v.ALT[0],
                "variant": v.REF + str(v.POS) + v.ALT[0],
                "amino_acid_consequence": "None",
                "nsp_aa_change": "None",
                "bcsq_nt_notation": "None",
                "bcsq_aa_notation": "None",
                "type_of_variant": v.var_type,
                "type_of_change": "None",
                "variant_status": variant_status,
                "persistence_status": persistent_status,
                "presence_absence": " / ".join(presence_absence),
                "first_appearance": first_appearance,
                "last_appearance": last_appearance,
                "overall_variant_qc": overall_variant_qc,
                "per_sample_variant_qc": " / ".join(depths_qc["variant_qc"]),
                "aa1_total_properties": "None",
                "aa2_total_properties": "None",
                "aa1_unique_properties": "None",
                "aa2_unique_properties": "None",
                "aa1_weight": "None",
                "aa2_weight": "None",
                "weight_difference": "None",
                "alt_freq": " / ".join(allele_freqs),
                "variant_depth": " / ".join(depths_qc["variant_depth"]),
                "variant_site_depth": " / ".join(depths_qc["site_depth"]),
                "variant_window_depth": " / ".join(depths_qc["window_depth"]),
                "samples": " / ".join(samples),
                "total_genome_coverage": " / ".join(total_cov_list),
            }
            results.append(result)
    results = pd.DataFrame(results)
    return results


def process_joint_variants(path):
    tab = pd.read_csv(path)
    tab = tab.assign(joint_variant=False)
    idx = tab[tab["bcsq_aa_notation"].str.startswith("@")].index
    for i in idx:
        # find the main pos and main index number
        main_pos = int(tab.loc[i]["bcsq_aa_notation"].replace("@", ""))
        j = tab[tab["start"] == main_pos].index[0]
        # update the joint variant key
        tab.at[i, "joint_variant"] = True
        tab.at[j, "joint_variant"] = True
        # assign the main position's value to the other variant
        tab.at[i, "gene"] = tab.at[j, "gene"]
        tab.at[i, "amino_acid_consequence"] = tab.at[j, "amino_acid_consequence"]
        tab.at[i, "nsp_aa_change"] = tab.at[j, "nsp_aa_change"]
        tab.at[i, "bcsq_nt_notation"] = tab.at[j, "bcsq_nt_notation"]
        tab.at[i, "bcsq_aa_notation"] = tab.at[j, "bcsq_aa_notation"]
        tab.at[i, "type_of_change"] = "joint_" + tab.at[j, "type_of_change"]
        tab.at[j, "type_of_change"] = "joint_" + tab.at[j, "type_of_change"]
        tab.at[i, "aa1_total_properties"] = tab.at[j, "aa1_total_properties"]
        tab.at[i, "aa2_total_properties"] = tab.at[j, "aa2_total_properties"]
        tab.at[i, "aa1_unique_properties"] = tab.at[j, "aa1_unique_properties"]
        tab.at[i, "aa2_unique_properties"] = tab.at[j, "aa2_unique_properties"]
        tab.at[i, "aa1_weight"] = tab.at[j, "aa1_weight"]
        tab.at[i, "aa2_weight"] = tab.at[j, "aa2_weight"]
        tab.at[i, "weight_difference"] = tab.at[j, "weight_difference"]
    tab.to_csv(path, index=None)
    return tab


def get_total_passage_mutations(row):
    colsum = sum(
        [
            int(x.replace("N", "0")) if x == "N" else int(x.replace("Y", "1"))
            for x in row
        ]
    )
    return colsum


def generate_cumulative_lineplot(table, pname, sample_number_list, outname):
    title = f"{pname}: cumulative mutations" if pname != "" else "Cumulative mutations"
    df = pd.DataFrame([x.split(" / ") for x in table["presence_absence"]]).transpose()
    plot_df = df
    plot_df["Cumulative Mutations (Total)"] = df.apply(
        get_total_passage_mutations, axis=1
    )
    mask = df.iloc[0] == "N"
    ndf = df.loc[:, mask]
    new_muts = ndf.apply(get_total_passage_mutations, axis=1)
    plot_df["Cumulative Mutations (New)"] = new_muts
    plot_df["Passage"] = sample_number_list
    plot_df = plot_df[
        ["Passage", "Cumulative Mutations (Total)", "Cumulative Mutations (New)"]
    ]
    plot_df = pd.melt(plot_df, ["Passage"], var_name="Type")
    f = sns.relplot(
        data=plot_df,
        x="Passage",
        y="value",
        hue="Type",
        kind="line",
        height=3.5,
        aspect=1.5,
    ).set(ylabel="Number of Mutations", xlabel="Passage", ylim=(0, None))
    f.fig.subplots_adjust(top=0.9)
    f.fig.suptitle(t=title, weight="bold")
    plt.savefig(outname, dpi=300)


def generate_gene_table(table):
    ref_gene_lengths = {
        "5' UTR": 265,
        "ORF1ab": 21291,
        "S": 3822,
        "ORF3a": 828,
        "E": 228,
        "M": 669,
        "ORF6": 186,
        "ORF7a": 366,
        "ORF7b": 132,
        "ORF8": 366,
        "N": 1260,
        "ORF10": 117,
        "3' UTR": 228,
        "INTERGENIC": 1,
    }
    genes = table["gene"].unique()
    change_types = table["type_of_change"].unique()
    result = []
    for gene in genes:
        if gene == "INTERGENIC":
            continue
        df = table[table["gene"] == gene]
        gene_result = []
        for ctype in change_types:
            if ctype == "None":
                continue
            out = {"gene": gene}
            out["type"] = ctype
            count = len(df[df["type_of_change"] == ctype])
            out["number"] = count
            gene_result.append(out)
        new_muts = {
            "gene": gene,
            "type": "new_mutations",
            "number": len([x for x in df["presence_absence"] if x.startswith("N")]),
        }
        persistent_muts = {
            "gene": gene,
            "type": "persistent_mutations",
            "number": len(
                [
                    x
                    for x in df["presence_absence"]
                    if x.startswith("N") and x.endswith("Y")
                ]
            ),
        }
        total = {
            "gene": gene,
            "type": "total",
            "number": sum([x.get("number") for x in gene_result]),
        }
        per_kb = {
            "gene": gene,
            "type": "mutations_per_kb",
            "number": (total.get("number") / ref_gene_lengths.get(gene)) * 1000,
        }
        new_per_kb = {
            "gene": gene,
            "type": "new_mutations_per_kb",
            "number": (new_muts.get("number") / ref_gene_lengths.get(gene)) * 1000,
        }
        gene_result.append(total)
        gene_result.append(new_muts)
        gene_result.append(persistent_muts)
        gene_result.append(per_kb)
        gene_result.append(new_per_kb)
        result.append(pd.DataFrame(gene_result))
    return pd.concat(result)


def plot_gene_table(gene_table, pname, outdir):
    g = sns.catplot(
        x="gene",
        y="number",
        col="type",
        col_wrap=3,
        data=gene_table,
        kind="bar",
        height=4,
        aspect=1.2,
    )
    g.set_axis_labels("", "Number of Mutations")
    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle(f"{pname}", weight="bold")
    for ax in g.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(45)
    for ax, title in zip(g.fig.axes, list(gene_table.type.unique())):
        ax.set_title(string.capwords(title.replace("_", " ")))
    plt.savefig(os.path.join(outdir, "mutations_per_gene.pdf"), dpi=300)


def search_pokay(table, pokay, outdir):
    print("Searching 'pokay' database to find significant mutations...\n")
    collector = []
    for _, row in table.iterrows():
        gene = row["gene"]
        mut = row["amino_acid_consequence"]
        gdf = pokay[(pokay.gene == gene) & (pokay.mutation.str.contains(mut))]
        if len(gdf) == 0:
            row["key_mutation"] = None
            row["database_mutation_string"] = None
            row["category"] = None
            row["prior_information"] = None
            row["reference"] = None
            collector.append(row)
        elif len(gdf) >= 1:
            for i in range(len(gdf)):
                r = row.copy()
                r["key_mutation"] = True
                r["database_mutation_string"] = gdf.iloc[i]["mutation"]
                r["category"] = gdf.iloc[i]["category"]
                r["prior_information"] = re.sub("^ +", "", gdf.iloc[i]["information"])
                r["reference"] = gdf.iloc[i]["reference"]
                collector.append(r)
    # make output with all info
    merged = pd.DataFrame(collector)
    key_muts = merged.loc[merged["key_mutation"] == True]
    to_keep = [
        "gene",
        "amino_acid_consequence",
        "database_mutation_string",
        "category",
        "prior_information",
        "reference",
    ]
    if "name" in key_muts.columns:
        reformatted1 = key_muts.groupby(to_keep, as_index=False)["name"].agg(", ".join)
        reformatted1.insert(0, "name", reformatted1.pop("name"))
    else:
        reformatted1 = key_muts[to_keep]
    fname1 = os.path.join(outdir, "pokay_database_hits.full.csv")
    reformatted1.to_csv(fname1, index=None)
    # make output with concise info
    merged = pd.DataFrame(collector)
    key_muts = merged.loc[merged["key_mutation"] == True]
    to_keep = ["gene", "amino_acid_consequence", "category", "reference", "name"]
    key_muts = key_muts[to_keep].drop_duplicates()
    to_keep = ["gene", "amino_acid_consequence", "category", "reference"]
    if "name" in key_muts.columns:
        reformatted = key_muts.groupby(to_keep, as_index=False)["name"].agg(", ".join)
        reformatted.insert(0, "name", reformatted.pop("name"))
        to_keep = ["name", "gene", "amino_acid_consequence", "category"]
    else:
        reformatted = key_muts[to_keep]
        to_keep = ["gene", "amino_acid_consequence", "category"]
    reformatted = reformatted.groupby(to_keep, as_index=False)["reference"].agg(
        " ; ".join
    )
    reformatted["reference"] = [
        re.sub(" ; ", "; ", x) for x in reformatted["reference"]
    ]
    reformatted["reference"] = [re.sub("\]", "", x) for x in reformatted["reference"]]
    l = []
    for x in reformatted["reference"]:
        y = x.split("; ")
        y = list(set(y))
        l.append("; ".join(y))
    reformatted["reference"] = l
    fname = os.path.join(outdir, "pokay_database_hits.concise.csv")
    reformatted.drop_duplicates().to_csv(fname, index=None)
    # print stats
    genes = len(reformatted["amino_acid_consequence"].unique())
    hits = len(reformatted["category"])
    print(f"{genes} new amino acid changes had hits to {hits} categories")
    print(f"See: {fname}")


def main(sysargs=sys.argv[1:]):
    print(logo)
    parser = argparse.ArgumentParser(
        description="vartracker: track the persistence (or not) or mutations during long-term passaging",
        usage="""python3 vartracker.py [options] <input_csv> """,
        formatter_class=FlexiFormatter,
        epilog="""
    The input CSV file should have four columns: 
    1. 'vcf': full paths to bgzipped vcf files containing the called variants for each sample of interest
    2. 'coverage': full paths to files containing the genome coverage for each sample. The coverage files should be in the format output by `bedtools genomecov` (e.g., bedtools genomecov -ibam input.bam -d > coverage.txt), whereby the columns are (in order): the reference name, the 1-based position in the genome, and the depth of coverage at that position. 
    3. 'sample_name': the name of the sample in the given VCF file
    4. 'sample_number': the sample number. In a longitudinal comparison like long-term passaging, the numbers in the column might go 0, 1, 2, ..., 15.
    
    Note: the order of values in the input CSV file matters, dictating the order of results in the output CSV file. 

    If 'new mutations' are to be searched against the functional 'pokay' database (default) (see vartracker repository README), the path to
    the parsed pokay csv file must be provided.
    """,
    )

    parser.add_argument("input_csv", nargs="*", help="Input CSV file. See below.")
    parser.add_argument(
        "-a",
        "--annotation",
        action="store",
        required=False,
        default=os.path.join(thisdir, "util", "NC_045512.gff3"),
        help="Annotations to use in GFF3 format",
    )
    parser.add_argument(
        "-c",
        "--chrom_map",
        action="store",
        required=False,
        default=os.path.join(thisdir, "util", "chrom_map.txt"),
        help="[DEPRECATED] File with map for renaming chrom in VCF files",
    )
    parser.add_argument(
        "-m",
        "--min_snv_freq",
        action="store",
        required=False,
        default=0.03,
        help="Minimum allele frequency of SNV variants to keep",
    )
    parser.add_argument(
        "-M",
        "--min_indel_freq",
        action="store",
        required=False,
        default=0.1,
        help="Minimum allele frequency of indel variants to keep",
    )
    parser.add_argument(
        "-n",
        "--name",
        action="store",
        required=False,
        default=None,
        help="Optional: add a column to results with the name specified here",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        required=False,
        default=thisdir,
        help="Output directory",
    )
    parser.add_argument(
        "-p",
        "--pokay",
        action="store",
        required=False,
        default=None,
        help="Path to csv file of the parsed 'pokay' database",
    )
    parser.add_argument(
        "-f",
        "--filename",
        action="store",
        required=False,
        default="results.csv",
        help="Output file name",
    )
    parser.add_argument(
        "-r",
        "--reference",
        action="store",
        required=False,
        default=os.path.join(thisdir, "util", "NC_045512.fasta"),
        help="Reference genome",
    )
    parser.add_argument(
        "--passage_cap",
        action="store",
        help="Cap the number of passages at this number",
        default=None,
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print commands being run for debugging",
        default=False,
    )
    parser.add_argument(
        "--keep_temp",
        action="store_true",
        help="Keep temporary files for debugging",
        default=False,
    )
    parser.add_argument(
        "--rename_chrom",
        action="store_true",
        help="[DEPRECATED] Rename chromosomes in VCF files. Useful if different reference names were given during VCF calling. Uses 'chrom_map' file. See bcftools/samtools usage guide.",
        default=False,
    )
    parser.add_argument(
        "--skip_search",
        action="store_true",
        help="Skip literature search for 'new mutations'.",
        default=False,
    )

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # read in input CSV
    try:
        input_file = pd.read_csv("".join(args.input_csv))
    except:
        parser.print_help()
        print("\n###########\n   ERROR\n###########")
        print(f"Could not read input file: {''.join(args.input_csv)}")
        print("Is the path correct?")
        print("Did you incorrectly specify multiple files?")
        print("See help above for more information")
        print("\nQuitting\n")
        sys.exit(1)

    if args.passage_cap is not None:
        input_file = input_file[input_file["sample_number"] <= int(args.passage_cap)]

    if not args.skip_search:
        if args.pokay is None:
            parser.print_help()
            print("\n###########\n   ERROR\n###########")
            print(
                "No path to 'pokay' CSV file provided yet the search is not being skipped"
            )
            print("See documentation")
            sys.exit(1)
        else:
            pokay = pd.read_csv(args.pokay)

    # set up variables
    vcfs = list(input_file["vcf"])
    sample_names = list(input_file["sample_name"])
    covs = list(input_file["coverage"])
    rename_chrom = True if args.rename_chrom else False
    reference = args.reference
    annotation = args.annotation
    chrom_map = args.chrom_map
    outdir = args.outdir
    debug = args.debug

    # make outdir if necessary
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # set up tsv outfile to store results
    outfile = os.path.join(outdir, args.filename)
    with open(outfile, "w") as output_handle:
        output_handle.write(
            "reference\tposition\tgene\tref_base\talt_base\tamino_acid_consequence\tbcsq_aa_notation\ttype_of_change\talt_freq\tpresence_absence\tsamples\n"
        )

    # make a tempdir for file processing
    tempdir = os.path.join(outdir, f"tempdir_{randrange(10,10000):04}")
    os.makedirs(tempdir)

    print("Pre-processing VCF files for compatibility...")
    # make sure they are all in a compatible format with the same chrom name etc.
    for vcf, sample in zip(vcfs, sample_names):
        format_vcf(
            vcf,
            sample,
            tempdir,
            rename_chrom,
            chrom_map,
            args.min_snv_freq,
            args.min_indel_freq,
            reference,
            annotation,
            debug,
        )

    new_vcfs = [
        re.sub(".vcf.gz", ".csq.vcf.gz", os.path.join(tempdir, os.path.basename(x)))
        for x in vcfs
    ]

    # write the paths to each new vcf in a file for use with  bcftools merge in the get_consequences function
    with open(os.path.join(tempdir, "vcf_list.txt"), "w") as output_handle:
        for vcf in new_vcfs:
            output_handle.write("%s\n" % vcf)

    # write the sample names to file for use with bcftools reheader in the merge_consequences function
    sample_names = os.path.join(tempdir, "sample_names.txt")
    with open(sample_names, "w") as output_handle:
        for sample_name in list(input_file["sample_name"]):
            output_handle.write("%s\n" % sample_name)

    # annotate the functional consequences using bcftools csq and merge VCFs
    print("Annotating results...")
    csq_file = os.path.join(tempdir, "vcf_annotated.vcf")
    merge_consequences(tempdir, csq_file, sample_names, debug)

    # summarise the called variants and add info longitudinally
    print("Summarising results...")
    table = process_vcf(csq_file, covs)

    # add in passage line/longitudinal name if desired
    if args.name is not None:
        table["name"] = args.name
        pname = args.name
    else:
        pname = ""

    # write to file
    table.to_csv(outfile, index=None)

    # process the joint variants from bcftools csq
    # deals with consequences starting with "@" for joint variants
    table = process_joint_variants(outfile)

    # plot cumulative results
    print("Plotting results...")
    generate_cumulative_lineplot(
        table,
        pname,
        input_file["sample_number"],
        os.path.join(outdir, "cumulative_mutations.pdf"),
    )
    gene_table = generate_gene_table(table)
    plot_gene_table(gene_table, pname, outdir)

    # write out specialised tables
    new_mutations = table[(table.presence_absence.str[0] != "Y")][
        ["gene", "variant", "amino_acid_consequence", "nsp_aa_change"]
    ].reset_index(drop=True)
    new_mutations.to_csv(os.path.join(outdir, "new_mutations.csv"), index=None)
    persistent_mutations = table[(table.persistence_status == "new_persistent")][
        ["gene", "variant", "amino_acid_consequence", "nsp_aa_change"]
    ].reset_index(drop=True)
    persistent_mutations.to_csv(
        os.path.join(outdir, "persistent_new_mutations.csv"), index=None
    )

    # search for new mutations in functional database
    if not args.skip_search:
        # search for importance of new mutations
        search_pokay(table, pokay, outdir)

    if not args.keep_temp:
        shutil.rmtree(tempdir)
    print(f"\nFinished: find results in {outdir}\n")


if __name__ == "__main__":
    main()
