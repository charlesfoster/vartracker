#!/usr/bin/env python3
"""Simulate longitudinal VCF + coverage data against the PAO1 reference bundle
for vartracker bacterial-scale validation (Scenario A ~80 variants, Scenario B ~1000).
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

N_TIMEPOINTS = 6
BASES = ["A", "C", "G", "T"]

PATTERNS = [
    "original_retained",
    "original_lost",
    "new_persistent",
    "new_transient",
    "original_intermittent",
    "new_intermittent",
]
PATTERN_WEIGHTS = [0.30, 0.15, 0.25, 0.15, 0.08, 0.07]
MIN_INTERMITTENT = 3


def parse_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:].split()[0]
            else:
                seq_chunks.append(line)
    return header, "".join(seq_chunks).upper()


def parse_cds(gff3_path, contig):
    cds = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != contig or parts[2] != "CDS":
                continue
            start, end, strand, attrs = int(parts[3]), int(parts[4]), parts[6], parts[8]
            gene = None
            for kv in attrs.split(";"):
                if kv.startswith("gene="):
                    gene = kv.split("=", 1)[1]
                    break
            if gene is None or end - start < 20:
                continue
            cds.append({"gene": gene, "start": start, "end": end, "strand": strand})
    return cds


def zipf_weights(n, rng, exponent=1.1):
    order = rng.permutation(n)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    w = 1.0 / (ranks ** exponent)
    return w / w.sum()


def choose_positions(cds_list, n_variants, ref_seq, rng, indel_fraction=0.08, biased=False):
    """Pick n_variants distinct (pos, ref, alt, gene, var_type) tuples from CDS regions."""
    genes = np.array([c["gene"] for c in cds_list])
    unique_genes, gene_index = np.unique(genes, return_inverse=True)
    # map each unique gene name -> list of cds_list indices (handles duplicate-name genes e.g. speA)
    gene_to_cds_idx = {g: [] for g in unique_genes}
    for i, g in enumerate(genes):
        gene_to_cds_idx[g].append(i)

    if biased:
        weights = zipf_weights(len(unique_genes), rng)
    else:
        weights = np.full(len(unique_genes), 1.0 / len(unique_genes))

    used_positions = set()
    variants = []
    attempts = 0
    max_attempts = n_variants * 50

    while len(variants) < n_variants and attempts < max_attempts:
        attempts += 1
        gene_choice = rng.choice(unique_genes, p=weights)
        cds_idx = rng.choice(gene_to_cds_idx[gene_choice])
        cds = cds_list[cds_idx]
        is_indel = rng.random() < indel_fraction
        span = rng.integers(1, 3) if is_indel else 0
        lo, hi = cds["start"] + 3, cds["end"] - 3 - span
        if hi <= lo:
            continue
        pos = int(rng.integers(lo, hi))
        span_positions = set(range(pos, pos + span + 2))
        if span_positions & used_positions:
            continue

        ref_base = ref_seq[pos - 1]
        if ref_base not in BASES:
            continue

        if not is_indel:
            alt_base = rng.choice([b for b in BASES if b != ref_base])
            ref, alt, var_type = ref_base, alt_base, "snp"
        else:
            if rng.random() < 0.5:
                # insertion
                inserted = "".join(rng.choice(BASES, size=span))
                ref, alt, var_type = ref_base, ref_base + inserted, "indel"
            else:
                # deletion: anchor base + (span) deleted ref bases
                del_seq = ref_seq[pos - 1: pos - 1 + span + 1]
                if len(del_seq) < span + 1 or any(b not in BASES for b in del_seq):
                    continue
                ref, alt, var_type = del_seq, del_seq[0], "indel"

        used_positions |= span_positions
        variants.append(
            {
                "chrom": None,  # filled by caller
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "gene": gene_choice,
                "var_type": var_type,
            }
        )

    if len(variants) < n_variants:
        print(
            f"WARNING: only placed {len(variants)}/{n_variants} variants after {attempts} attempts",
            file=sys.stderr,
        )
    return variants


def gen_presence_af(pattern, var_type, rng):
    af_floor = 0.12 if var_type == "indel" else 0.05
    af_ceiling = 0.95
    present = [False] * N_TIMEPOINTS

    if pattern == "original_retained":
        present = [True] * N_TIMEPOINTS
        base = rng.uniform(0.2, 0.5)
        af_vals = [np.clip(base + rng.uniform(-0.05, 0.05), af_floor, af_ceiling) for _ in range(N_TIMEPOINTS)]

    elif pattern == "original_lost":
        last_present = int(rng.integers(1, 5))  # 1..4 inclusive -> absent by tp5
        present = [i <= last_present for i in range(N_TIMEPOINTS)]
        n_present = last_present + 1
        af_vals = list(np.clip(np.linspace(0.55, 0.15, n_present) + rng.uniform(-0.03, 0.03, n_present), af_floor, af_ceiling))

    elif pattern == "new_persistent":
        first_present = int(rng.integers(1, 5))  # 1..4
        present = [i >= first_present for i in range(N_TIMEPOINTS)]
        n_present = N_TIMEPOINTS - first_present
        af_vals = list(np.clip(np.linspace(0.12, 0.85, n_present) + rng.uniform(-0.03, 0.03, n_present), af_floor, af_ceiling))

    elif pattern == "new_transient":
        first_present = int(rng.integers(1, 4))  # 1..3
        last_present = int(rng.integers(first_present, min(first_present + 3, 5)))  # < 5
        present = [first_present <= i <= last_present for i in range(N_TIMEPOINTS)]
        n_present = last_present - first_present + 1
        base = rng.uniform(0.15, 0.4)
        af_vals = list(np.clip([base + rng.uniform(-0.05, 0.05) for _ in range(n_present)], af_floor, af_ceiling))

    elif pattern == "original_intermittent":
        gap_start = int(rng.integers(1, 4))  # 1..3
        gap_len = int(rng.integers(1, 5 - gap_start))  # keep gap inside 1..4, tp5 stays present
        present = [True] * N_TIMEPOINTS
        for i in range(gap_start, min(gap_start + gap_len, 5)):
            present[i] = False
        present[0] = True
        present[5] = True
        n_present = sum(present)
        base = rng.uniform(0.2, 0.45)
        af_vals = list(np.clip([base + rng.uniform(-0.05, 0.05) for _ in range(n_present)], af_floor, af_ceiling))

    elif pattern == "new_intermittent":
        first_present = int(rng.integers(1, 3))  # 1..2
        gap_start = int(rng.integers(first_present + 1, 5))  # inside 1..4
        gap_len = int(rng.integers(1, max(2, 5 - gap_start)))
        present = [False] * N_TIMEPOINTS
        for i in range(first_present, N_TIMEPOINTS):
            present[i] = True
        for i in range(gap_start, min(gap_start + gap_len, 5)):
            present[i] = False
        present[0] = False
        present[5] = True
        n_present = sum(present)
        base = rng.uniform(0.15, 0.4)
        af_vals = list(np.clip([base + rng.uniform(-0.05, 0.05) for _ in range(n_present)], af_floor, af_ceiling))

    else:
        raise ValueError(pattern)

    # sanity checks on invariants
    if pattern.startswith("original"):
        assert present[0] is True
    else:
        assert present[0] is False
    if pattern in ("original_retained", "new_persistent", "original_intermittent", "new_intermittent"):
        assert present[5] is True
    else:
        assert present[5] is False

    af_iter = iter(af_vals)
    af_full = [round(float(next(af_iter)), 4) if p else None for p in present]
    dp_full = [int(rng.integers(80, 200)) if p else None for p in present]
    return present, af_full, dp_full


def assign_patterns(n, rng):
    patterns = list(rng.choice(PATTERNS, size=n, p=PATTERN_WEIGHTS))
    for target in ("original_intermittent", "new_intermittent"):
        count = patterns.count(target)
        if count < MIN_INTERMITTENT:
            donor_pool = [i for i, p in enumerate(patterns) if p not in ("original_intermittent", "new_intermittent")]
            n_needed = MIN_INTERMITTENT - count
            idx_to_convert = rng.choice(donor_pool, size=n_needed, replace=False)
            for idx in idx_to_convert:
                patterns[idx] = target
    return patterns


def simulate_scenario(name, n_variants, cds_list, contig, ref_seq, seed, biased, outdir: Path):
    rng = np.random.default_rng(seed)
    variants = choose_positions(cds_list, n_variants, ref_seq, rng, biased=biased)
    for v in variants:
        v["chrom"] = contig

    patterns = assign_patterns(len(variants), rng)

    records = []
    for v, pattern in zip(variants, patterns):
        present, af, dp = gen_presence_af(pattern, v["var_type"], rng)
        variant_status = "original" if pattern.startswith("original") else "new"
        records.append(
            {
                **v,
                "pattern": pattern,
                "variant_status": variant_status,
                "present": present,
                "af": af,
                "dp": dp,
            }
        )

    records.sort(key=lambda r: r["pos"])

    # sanity: no duplicate positions
    positions = [r["pos"] for r in records]
    assert len(positions) == len(set(positions)), "duplicate positions generated"

    outdir.mkdir(parents=True, exist_ok=True)
    vcf_dir = outdir / "vcfs"
    vcf_dir.mkdir(exist_ok=True)

    for tp in range(N_TIMEPOINTS):
        rows = [r for r in records if r["present"][tp]]
        rows.sort(key=lambda r: r["pos"])
        raw_path = vcf_dir / f"pao1_tp{tp}.raw.vcf"
        with open(raw_path, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            fh.write(f"##contig=<ID={contig}>\n")
            fh.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw Depth">\n')
            fh.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n')
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for r in rows:
                af_val = r["af"][tp]
                dp_val = r["dp"][tp]
                qual = 60
                fh.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t{qual}\tPASS\tDP={dp_val};AF={af_val}\n")

        final_gz = vcf_dir / f"pao1_tp{tp}.vcf.gz"
        sort_cmd = f"bcftools sort -Ov {raw_path} 2>/dev/null | bgzip > {final_gz}"
        subprocess.run(sort_cmd, shell=True, check=True)
        subprocess.run(["bcftools", "index", "-f", str(final_gz)], check=True)
        raw_path.unlink()

    # ground truth table for cross-checking step 4
    gt_rows = []
    for r in records:
        gt_rows.append(
            {
                "chrom": r["chrom"],
                "pos": r["pos"],
                "ref": r["ref"],
                "alt": r["alt"],
                "gene": r["gene"],
                "var_type": r["var_type"],
                "variant_status": r["variant_status"],
                "intended_pattern": r["pattern"],
                "presence": "/".join("Y" if p else "N" for p in r["present"]),
                "af": "/".join("" if a is None else str(a) for a in r["af"]),
            }
        )
    gt_df = pd.DataFrame(gt_rows).sort_values(["gene", "pos"])
    gt_df.to_csv(outdir / f"{name}_ground_truth.csv", index=False)

    # input CSV for vartracker vcf mode
    csv_rows = []
    for tp in range(N_TIMEPOINTS):
        csv_rows.append(
            {
                "sample_name": f"pao1_tp{tp}",
                "sample_number": tp,
                "reads1": "",
                "reads2": "",
                "bam": "",
                "vcf": str((vcf_dir / f"pao1_tp{tp}.vcf.gz").resolve()),
                "coverage": str((outdir.parent / "coverage" / f"pao1_tp{tp}_depth.txt").resolve()),
            }
        )
    pd.DataFrame(csv_rows).to_csv(outdir / f"{name}_input.csv", index=False)

    # summary counts
    pattern_counts = pd.Series([r["pattern"] for r in records]).value_counts()
    gene_counts = pd.Series([r["gene"] for r in records]).value_counts()
    print(f"[{name}] {len(records)} distinct variants placed across {gene_counts.shape[0]} genes")
    print(f"[{name}] pattern mix:\n{pattern_counts}")
    print(f"[{name}] top genes by variant count:\n{gene_counts.head(10)}")
    return records


def generate_coverage(contig, genome_length, coverage_dir: Path, seed=123):
    coverage_dir.mkdir(parents=True, exist_ok=True)
    pos = np.arange(1, genome_length + 1, dtype=np.int64)
    for tp in range(N_TIMEPOINTS):
        out_path = coverage_dir / f"pao1_tp{tp}_depth.txt"
        if out_path.exists():
            continue
        rng = np.random.default_rng(seed + tp)
        depth = rng.integers(80, 151, size=genome_length, dtype=np.int32)
        df = pd.DataFrame({"chrom": contig, "pos": pos, "depth": depth})
        df.to_csv(out_path, sep="\t", header=False, index=False)
        print(f"wrote coverage file {out_path} ({genome_length} lines)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--base-outdir", required=True)
    args = ap.parse_args()

    base = Path(args.base_outdir)
    contig, ref_seq = parse_fasta(args.fasta)
    print(f"Loaded reference contig {contig}, length {len(ref_seq)}")
    cds_list = parse_cds(args.gff3, contig)
    print(f"Parsed {len(cds_list)} CDS features on {contig}")

    generate_coverage(contig, len(ref_seq), base / "coverage")

    simulate_scenario(
        "scenario_a", 80, cds_list, contig, ref_seq, seed=42, biased=False, outdir=base / "scenario_a"
    )
    simulate_scenario(
        "scenario_b", 1000, cds_list, contig, ref_seq, seed=43, biased=True, outdir=base / "scenario_b"
    )


if __name__ == "__main__":
    main()
