"""
Analysis and plotting functionality for vartracker.

Contains functions for generating plots and analyzing mutation patterns.
"""

import os
import string
from typing import Dict, List, Optional, Sequence, Tuple, Union, cast

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from .constants import REF_GENE_LENGTHS, NSP_LENGTHS, NSPS

# Global configuration for plotting
plt.rcdefaults()
mpl.rcParams["pdf.fonttype"] = 42


def process_joint_variants(path):
    """
    Process joint variants from bcftools csq output.

    Args:
        path (str): Path to CSV file with variants

    Returns:
        pd.DataFrame: Updated DataFrame with joint variant information
    """
    try:
        tab = pd.read_csv(path)
    except Exception as e:
        raise RuntimeError(f"Error reading variants file {path}: {str(e)}")

    tab = tab.assign(joint_variant=False)

    # Get indices where bcsq_aa_notation starts with "@"
    idx = tab[tab["bcsq_aa_notation"].str.startswith("@", na=False)].index

    # If no joint variants, return the table as is
    if idx.empty:
        return tab

    for i in idx:
        try:
            # Find the main pos and main index number
            main_pos = int(tab.loc[i]["bcsq_aa_notation"].replace("@", ""))
            j = tab[tab["start"] == main_pos].index[0]

            # Update the joint variant key
            tab.at[i, "joint_variant"] = True
            tab.at[j, "joint_variant"] = True

            # Assign the main position's value to the other variant
            for col in [
                "gene",
                "amino_acid_consequence",
                "nsp_aa_change",
                "bcsq_nt_notation",
                "bcsq_aa_notation",
                "aa1_total_properties",
                "aa2_total_properties",
                "aa1_unique_properties",
                "aa2_unique_properties",
                "aa1_weight",
                "aa2_weight",
                "weight_difference",
            ]:
                tab.at[i, col] = tab.at[j, col]

            # Update type of change for joint variants
            tab.at[i, "type_of_change"] = "joint_" + tab.at[j, "type_of_change"]
            tab.at[j, "type_of_change"] = "joint_" + tab.at[j, "type_of_change"]

        except (ValueError, IndexError, KeyError) as e:
            print(f"Warning: Could not process joint variant at index {i}: {str(e)}")
            continue

    tab.to_csv(path, index=None)
    return tab


def get_total_passage_mutations(row):
    """
    Calculate total number of mutations present across passages.

    Args:
        row: Pandas Series with presence/absence data

    Returns:
        int: Total count of mutations
    """
    return sum([int(x.replace("N", "0").replace("Y", "1")) for x in row])


def generate_cumulative_lineplot(table, pname, sample_number_list, outname):
    """
    Generate cumulative mutations plot.

    Args:
        table (pd.DataFrame): Variants table
        pname (str): Project name for plot title
        sample_number_list: List of sample numbers
        outname (str): Output file path
    """
    try:
        title = (
            f"{pname}: cumulative mutations" if pname != "" else "Cumulative mutations"
        )
        df = pd.DataFrame(
            [x.split(" / ") for x in table["presence_absence"]]
        ).transpose()
        plot_df = df.copy()

        plot_df["Cumulative Mutations (Total)"] = df.apply(
            get_total_passage_mutations, axis=1
        )

        # Calculate new mutations (those not present in first sample)
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
        ).set(ylabel="Number of Mutations", xlabel="Sample", ylim=(0, None))

        ax = f.ax
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

        f.fig.subplots_adjust(top=0.9)
        f.fig.suptitle(t=title, weight="bold")
        plt.savefig(outname, dpi=300, bbox_inches="tight")
        plt.close()

    except Exception as e:
        raise RuntimeError(f"Error generating cumulative line plot: {str(e)}")


def generate_gene_table(
    table: pd.DataFrame, gene_lengths: Dict[str, int] | None = None
):
    """
    Generate gene-wise mutation statistics table.

    Args:
        table (pd.DataFrame): Variants table

    Returns:
        pd.DataFrame: Gene statistics table
    """
    # Combine gene sets: canonical genes + nsps
    source_gene_lengths = gene_lengths or REF_GENE_LENGTHS
    all_genes = list(source_gene_lengths.keys())
    include_nsps = gene_lengths is None
    nsp_products: List[str] = []
    if include_nsps:
        nsp_products = list(cast(Sequence[str], NSPS.get("product", [])))
        all_genes.extend(nsp_products)

    change_types = table["type_of_change"].unique()

    def make_gene_rows(gene, df):
        """Make rows for a single gene or NSP."""
        gene_result = []
        for ctype in change_types:
            if ctype == "None" or pd.isna(ctype):
                continue
            count = len(df[df["type_of_change"] == ctype])
            gene_result.append({"gene": gene, "type": ctype, "number": count})

        new_muts_count = (
            len([x for x in df["presence_absence"] if x.startswith("N")])
            if not df.empty
            else 0
        )

        persistent_count = (
            len(
                [
                    x
                    for x in df["presence_absence"]
                    if x.startswith("N") and x.endswith("Y")
                ]
            )
            if not df.empty
            else 0
        )

        total_count = sum(r["number"] for r in gene_result)

        length = source_gene_lengths.get(gene, NSP_LENGTHS.get(gene, 1))
        per_kb = (total_count / length) * 1000 if length else 0
        new_per_kb = (new_muts_count / length) * 1000 if length else 0

        gene_result.extend(
            [
                {"gene": gene, "type": "total", "number": total_count},
                {"gene": gene, "type": "new_mutations", "number": new_muts_count},
                {
                    "gene": gene,
                    "type": "persistent_mutations",
                    "number": persistent_count,
                },
                {"gene": gene, "type": "mutations_per_kb", "number": per_kb},
                {"gene": gene, "type": "new_mutations_per_kb", "number": new_per_kb},
            ]
        )

        return pd.DataFrame(gene_result)

    # Process original genes
    result = []
    genes_in_table = table["gene"].unique()

    for gene in genes_in_table:
        if gene == "INTERGENIC":
            continue
        df = table[table["gene"] == gene]
        result.append(make_gene_rows(gene, df))

    if include_nsps:
        table = table.copy()

        def assign_nsp(row):
            if (
                row["gene"] == "ORF1ab"
                and row["nsp_aa_change"]
                and ":" in str(row["nsp_aa_change"])
            ):
                return row["nsp_aa_change"].split(":")[0]
            return None

        table["nsp_gene"] = table.apply(assign_nsp, axis=1)
        orf_df = table[table["gene"] == "ORF1ab"]

        for nsp_gene_name in nsp_products:
            nsp_df = orf_df[orf_df["nsp_gene"] == nsp_gene_name]
            if not nsp_df.empty:
                result.append(make_gene_rows(nsp_gene_name, nsp_df))

    if not result:
        return pd.DataFrame()

    gene_table = pd.concat(result, ignore_index=True)

    # Create scaffold to ensure all genes and nsps appear even if zero variants
    all_types = gene_table["type"].unique()
    scaffold = pd.DataFrame(
        [(g, t) for g in all_genes for t in all_types], columns=["gene", "type"]
    )

    gene_table = scaffold.merge(gene_table, on=["gene", "type"], how="left")
    gene_table["number"] = gene_table["number"].fillna(0)

    return gene_table


def plot_gene_table(gene_table, pname, outdir):
    """
    Plot gene-wise mutation statistics.

    Args:
        gene_table (pd.DataFrame): Gene statistics table
        pname (str): Project name for plot title
        outdir (str): Output directory
    """
    try:
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

        # Rotate labels by 90 degrees
        for ax in g.axes.flat:
            for label in ax.get_xticklabels():
                label.set_rotation(90)

        # Make subplot titles nicer
        for ax, title in zip(g.fig.axes, list(gene_table["type"].unique())):
            # Convert title to string to handle numeric values and NaN
            if pd.isna(title):
                title_str = "None"
            else:
                title_str = str(title) if not isinstance(title, str) else title
            ax.set_title(string.capwords(title_str.replace("_", " ")))

        plt.savefig(
            os.path.join(outdir, "mutations_per_gene.pdf"), dpi=300, bbox_inches="tight"
        )
        plt.close()

    except Exception as e:
        raise RuntimeError(f"Error plotting gene table: {str(e)}")


def _build_gene_order_map(
    gene_lengths: Dict[str, int] | None = None,
    include_nsps: bool = True,
) -> Tuple[Dict[str, int], List[str]]:
    """Build an ordered mapping of genes (and NSPs) following 5' to 3'."""

    base_order = (
        list(gene_lengths.keys())
        if gene_lengths is not None
        else list(REF_GENE_LENGTHS.keys())
    )

    nsps = list(cast(Sequence[str], NSPS.get("product", [])))

    order: List[str] = []
    if "5' UTR" in base_order or gene_lengths is None:
        order.append("5' UTR")

    if include_nsps:
        order.extend(nsps)

    for gene in base_order:
        if gene in {"5' UTR", "ORF1ab"}:
            continue
        if gene not in order:
            order.append(gene)

    if include_nsps and "ORF1ab" in base_order:
        order.append("ORF1ab")

    seen = set()
    ordered = []
    for gene in order:
        if gene and gene not in seen:
            ordered.append(gene)
            seen.add(gene)

    gene_order_map = {gene: idx for idx, gene in enumerate(ordered)}
    if include_nsps and "nsp1" in gene_order_map:
        gene_order_map.setdefault("ORF1ab", gene_order_map["nsp1"])

    gene_order_map.setdefault("INTERGENIC", len(gene_order_map))

    return gene_order_map, ordered


def _map_orf1ab_position_to_nsp(aa_position: int) -> str:
    """Map an ORF1ab amino acid position to its NSP name."""
    nsp_products = cast(Sequence[str], NSPS.get("product", []))
    aa_starts = cast(Sequence[int], NSPS.get("aa_start", []))

    for name, start in zip(nsp_products, aa_starts):
        length = NSP_LENGTHS.get(str(name))
        if length is None:
            continue
        end = start + length - 1
        if start <= aa_position <= end:
            return name

    return "ORF1ab"


def _extract_numeric_position(value: str) -> Optional[int]:
    """Extract the numeric portion from a mutation string, if present."""
    if not isinstance(value, str):
        return None
    digits = "".join(ch for ch in value if ch.isdigit())
    return int(digits) if digits else None


def _extract_amino_acid_letter(value: str) -> str:
    """Return the first alphabetic character from a mutation string."""
    if not isinstance(value, str):
        return ""
    for ch in value:
        if ch.isalpha():
            return ch
    return ""


def _coerce_frequency(value: str) -> float:
    """Convert allele frequency strings to floats handling sentinel values."""
    if value is None:
        return 0.0
    if isinstance(value, float):
        return value
    text = str(value).strip()
    if text in {"", ".", "NA", "nan", "NaN"}:
        return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0


def _resolve_variant_labels(row) -> Tuple[str, str, str]:
    """Determine gene label and variant labels for heatmap plotting."""
    gene = getattr(row, "gene", "")
    amino_change = getattr(row, "amino_acid_consequence", "")
    nsp_change = getattr(row, "nsp_aa_change", "")
    change_type = str(getattr(row, "type_of_change", ""))

    gene_label = gene
    aa_label = amino_change

    if gene == "ORF1ab":
        if isinstance(nsp_change, str) and nsp_change and nsp_change != "ERROR":
            gene_part, _, aa_part = nsp_change.partition(":")
            if gene_part:
                gene_label = gene_part
            if aa_part:
                aa_label = aa_part
        else:
            position = _extract_numeric_position(amino_change)
            if position is not None:
                gene_label = _map_orf1ab_position_to_nsp(position)

    if not aa_label or str(aa_label) in {"", "None"}:
        if isinstance(nsp_change, str) and ":" in nsp_change:
            aa_label = nsp_change.split(":", 1)[1]
        elif amino_change:
            aa_label = amino_change
        else:
            aa_label = getattr(row, "variant", "")

    if change_type.lower() == "synonymous":
        position = _extract_numeric_position(aa_label)
        aa_letter = _extract_amino_acid_letter(aa_label)
        if position is not None and aa_letter:
            aa_label = f"{aa_letter}{position}="
        elif position is not None:
            aa_label = f"{position}="
        elif aa_label:
            aa_label = f"{aa_label}="

    aa_label = str(aa_label)
    gene_label = str(gene_label) if gene_label else gene

    trunc_limit = 10
    if len(aa_label) > trunc_limit * 2:
        head = aa_label[:trunc_limit]
        tail = aa_label[-trunc_limit:]
        middle = len(aa_label) - (trunc_limit * 2)
        if middle > 0:
            aa_label = f"{head}+{middle}{tail}"

    if aa_label.startswith(f"{gene_label}:"):
        base_label = aa_label
    elif aa_label:
        base_label = f"{gene_label}:{aa_label}"
    else:
        base_label = gene_label

    nuc_change = str(getattr(row, "variant", "")).strip()
    if nuc_change and nuc_change not in {"", "None"}:
        display_label = f"{base_label}\n({nuc_change})"
    else:
        display_label = base_label

    return gene_label, display_label, base_label


def _prepare_variant_heatmap_matrix(
    table: pd.DataFrame,
    sample_names: Sequence[str],
    min_snv_freq: float,
    min_indel_freq: float,
    gene_lengths: Dict[str, int] | None = None,
) -> pd.DataFrame:
    """Prepare a matrix of allele frequencies for heatmap plotting."""
    if not isinstance(table, pd.DataFrame) or table.empty:
        return pd.DataFrame(columns=list(sample_names))

    if "alt_freq" not in table.columns or "samples" not in table.columns:
        raise ValueError("table must contain 'alt_freq' and 'samples' columns")

    ordered_samples = [str(name) for name in sample_names]
    if not ordered_samples:
        ordered_samples = [
            s.strip() for s in str(table.iloc[0]["samples"]).split(" / ")
        ]

    use_nsps = True
    if gene_lengths is not None:
        use_nsps = False
        if "gene" in table.columns:
            genes_series = table["gene"].astype(str)
            use_nsps = (
                genes_series.isin(cast(Sequence[str], NSPS.get("product", []))).any()
                or genes_series.eq("ORF1ab").any()
            )

        if not use_nsps and "nsp_aa_change" in table.columns:
            use_nsps = table["nsp_aa_change"].astype(str).str.contains(":").any()

    gene_order_map, _ = _build_gene_order_map(gene_lengths, include_nsps=use_nsps)

    records: List[Dict[str, Union[str, float, int]]] = []
    seen_labels = set()

    for row in table.itertuples(index=False):
        gene_value = getattr(row, "gene", "")
        if str(gene_value) in {"5' UTR", "3' UTR", "INTERGENIC"}:
            continue

        gene_label, display_label, base_label = _resolve_variant_labels(row)

        if base_label in seen_labels:
            continue

        variant_type = str(getattr(row, "type_of_variant", "")).lower()

        sample_tokens = [
            s.strip() for s in str(getattr(row, "samples", "")).split(" / ")
        ]
        freq_tokens = [
            _coerce_frequency(token)
            for token in str(getattr(row, "alt_freq", "")).split(" / ")
        ]

        if not freq_tokens:
            continue

        max_freq = max(freq_tokens)

        threshold = min_snv_freq
        if "indel" in variant_type:
            threshold = min_indel_freq

        if max_freq < threshold:
            continue

        sample_freq_map = {
            sample: freq for sample, freq in zip(sample_tokens, freq_tokens)
        }
        row_values = [sample_freq_map.get(sample, 0.0) for sample in ordered_samples]

        record: Dict[str, Union[str, float, int]] = {
            "label": display_label,
            "base_label": base_label,
            "gene_order": gene_order_map.get(gene_label, len(gene_order_map)),
            "start": getattr(row, "start", 0),
        }
        for sample, value in zip(ordered_samples, row_values):
            record[sample] = value

        records.append(record)
        seen_labels.add(base_label)

    if not records:
        return pd.DataFrame(columns=ordered_samples)

    matrix = pd.DataFrame(records)
    matrix = matrix.sort_values(["gene_order", "start", "base_label"])
    matrix = matrix.drop(columns=["gene_order", "start"], errors="ignore")
    matrix = matrix.drop_duplicates(subset=["base_label"]).set_index("label")
    matrix = matrix.drop(columns=["base_label"], errors="ignore")

    matrix = matrix.loc[matrix.max(axis=1) > 0]

    return matrix.reindex(columns=ordered_samples).fillna(0.0)


def generate_variant_heatmap(
    table: pd.DataFrame,
    sample_names: Sequence[str],
    sample_numbers: Union[Sequence[int], Sequence[str]],
    outdir: str,
    project_name: str,
    min_snv_freq: float,
    min_indel_freq: float,
    gene_lengths: Dict[str, int] | None = None,
):
    """Generate a heatmap of variant allele frequencies across passages."""

    try:
        heatmap_data = _prepare_variant_heatmap_matrix(
            table, sample_names, min_snv_freq, min_indel_freq, gene_lengths
        )
        if heatmap_data.empty:
            print("No variant data available for heatmap; skipping plot.")
            return

        heatmap_path = os.path.join(outdir, "variant_allele_frequency_heatmap.pdf")

        fig_width = max(4, 1.2 * len(heatmap_data.columns))
        fig_height = max(4, 0.4 * len(heatmap_data))
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        plot_data = heatmap_data.replace(0, np.nan)
        cmap = sns.color_palette("viridis", as_cmap=True).copy()
        cmap.set_bad(color="#d9d9d9")
        sns.heatmap(
            plot_data,
            cmap=cmap,
            vmin=0,
            vmax=1,
            cbar_kws={"label": "Variant allele frequency"},
            linewidths=0.2,
            linecolor="white",
            ax=ax,
        )

        heatmap_title = (
            f"{project_name}: variant allele frequencies"
            if project_name
            else "Variant allele frequencies"
        )
        ax.set_title(heatmap_title, weight="bold")

        tick_labels = list(sample_names)
        ax.set_xticklabels(tick_labels, rotation=45, ha="right")
        ax.set_xlabel("Sample", fontweight="bold")
        ax.set_ylabel("Variant (5' → 3')", fontweight="bold")

        fig.tight_layout()
        fig.savefig(heatmap_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

    except Exception as exc:
        raise RuntimeError(f"Error generating variant heatmap: {exc}") from exc


def search_pokay(table, pokay, outdir, pokay_name, debug=False):
    """
    Search mutations against the pokay functional database.

    Args:
        table (pd.DataFrame): Variants table
        pokay (pd.DataFrame): Pokay database
        outdir (str): Output directory
        pokay_name (str): Name for output files
        debug (bool): Whether to print debug information

    Returns:
        pd.DataFrame: Full search results
    """
    print("Searching 'pokay' database to find significant mutations...\n")

    if debug:
        print(f"Input table shape: {table.shape}")
        print(f"Input table columns: {list(table.columns)}")
        print(f"Pokay database shape: {pokay.shape}")
        print(f"Pokay database columns: {list(pokay.columns)}")

    collector = []

    # Normalize pokay lookup columns to avoid nullable boolean masks
    pokay_gene = pokay.get("gene")
    pokay_mutation = pokay.get("mutation")

    if pokay_gene is None or pokay_mutation is None:
        raise ValueError("Pokay database must contain 'gene' and 'mutation' columns")

    pokay_gene_series = pokay_gene.astype("string").fillna("")
    pokay_mutation_series = pokay_mutation.astype("string").fillna("")

    for _, row in table.iterrows():
        gene = row["gene"]
        mut = row["amino_acid_consequence"]

        if gene == "ORF1ab" and ":" in str(row.get("nsp_aa_change", "")):
            search_gene = row["nsp_aa_change"].split(":")[0].split("_")[0]
            search_mut = row["nsp_aa_change"].split(":")[1]
        else:
            search_gene = gene
            search_mut = mut

        if debug:
            print(f"Searching: {search_gene} {search_mut}")

        try:
            search_gene_str = str(search_gene) if search_gene is not None else ""
            search_mut_str = str(search_mut) if search_mut is not None else ""

            gene_mask = (pokay_gene_series == search_gene_str).fillna(False)
            mutation_mask = pokay_mutation_series.str.contains(
                search_mut_str, regex=False, na=False
            ).fillna(False)

            combined_mask = pd.Series(gene_mask & mutation_mask, index=pokay.index)
            gdf = pokay.loc[combined_mask.astype(bool)]
        except Exception as e:
            if debug:
                print(f"Error searching pokay for {search_gene} {search_mut}: {str(e)}")
            gdf = pd.DataFrame()

        if len(gdf) == 0:
            row_copy = row.copy()
            row_copy["key_mutation"] = False
            row_copy["database_mutation_string"] = None
            row_copy["category"] = None
            row_copy["prior_information"] = None
            row_copy["reference"] = None
            collector.append(row_copy)
        else:
            for i in range(len(gdf)):
                row_copy = row.copy()
                row_copy["key_mutation"] = True
                row_copy["database_mutation_string"] = gdf.iloc[i]["mutation"]
                row_copy["category"] = gdf.iloc[i]["category"]
                row_copy["prior_information"] = str(gdf.iloc[i]["information"]).lstrip()
                row_copy["reference"] = gdf.iloc[i]["reference"]
                if gene == "ORF1ab":
                    row_copy["gene"] = search_gene
                    row_copy["amino_acid_consequence"] = search_mut
                collector.append(row_copy)

    # Create output with all info
    if debug:
        print(f"Collector length: {len(collector)}")
        if len(collector) > 0:
            print(
                f"First collector item keys: {list(collector[0].keys()) if collector else 'No items'}"
            )
        else:
            print("Collector is empty!")

    # Handle empty collector gracefully
    if len(collector) == 0:
        print("No mutations found to search against pokay database.")
        # Return empty DataFrames with expected structure
        empty_df = pd.DataFrame(
            columns=[
                "gene",
                "amino_acid_consequence",
                "database_mutation_string",
                "category",
                "prior_information",
                "reference",
            ]
        )
        fname1 = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.full.csv")
        empty_df.to_csv(fname1, index=None)
        fname = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.concise.csv")
        empty_df.drop(columns=["database_mutation_string", "prior_information"]).to_csv(
            fname, index=None
        )
        print("0 new amino acid changes had hits to 0 categories")
        print(f"See: {fname}")
        return empty_df

    merged = pd.DataFrame(collector).drop_duplicates()

    if debug:
        print(f"Merged DataFrame shape: {merged.shape}")
        print(f"Merged DataFrame columns: {list(merged.columns)}")

    # Check if key_mutation column exists
    if "key_mutation" not in merged.columns:
        print("Warning: 'key_mutation' column not found in merged DataFrame.")
        # Create empty results
        empty_df = pd.DataFrame(
            columns=[
                "gene",
                "amino_acid_consequence",
                "database_mutation_string",
                "category",
                "prior_information",
                "reference",
            ]
        )
        fname1 = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.full.csv")
        empty_df.to_csv(fname1, index=None)
        fname = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.concise.csv")
        empty_df.drop(columns=["database_mutation_string", "prior_information"]).to_csv(
            fname, index=None
        )
        print("0 new amino acid changes had hits to 0 categories")
        print(f"See: {fname}")
        return empty_df

    key_muts = merged.loc[merged["key_mutation"]]

    if key_muts.empty:
        empty_df = pd.DataFrame(
            columns=[
                "gene",
                "amino_acid_consequence",
                "database_mutation_string",
                "category",
                "prior_information",
                "reference",
            ]
        )
        fname1 = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.full.csv")
        empty_df.to_csv(fname1, index=None)
        fname = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.concise.csv")
        empty_df.drop(columns=["database_mutation_string", "prior_information"]).to_csv(
            fname, index=None
        )
        print("0 new amino acid changes had hits to 0 categories")
        print(f"See: {fname}")
        return empty_df

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

    fname1 = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.full.csv")
    reformatted1.to_csv(fname1, index=None)

    # Create concise output
    to_keep_concise = ["gene", "amino_acid_consequence", "category", "reference"]
    if "name" in key_muts.columns:
        to_keep_concise.insert(0, "name")

    key_muts_subset = key_muts[to_keep_concise].drop_duplicates()

    if "name" in key_muts_subset.columns:
        # First, group by everything except 'name' and aggregate 'name'
        temp_group_cols = [col for col in to_keep_concise if col != "reference"]
        reformatted = key_muts_subset.groupby(temp_group_cols, as_index=False).agg(
            {"name": ", ".join, "reference": " ; ".join}
        )
        reformatted.insert(0, "name", reformatted.pop("name"))
        group_cols = ["name", "gene", "amino_acid_consequence", "category"]
    else:
        reformatted = key_muts_subset
        group_cols = ["gene", "amino_acid_consequence", "category"]

    # Only do the second groupby if we didn't already aggregate references above
    if "name" not in key_muts_subset.columns:
        reformatted = reformatted.groupby(group_cols, as_index=False)["reference"].agg(
            " ; ".join
        )

    # Clean up references
    reformatted["reference"] = reformatted["reference"].str.replace(" ; ", "; ")
    reformatted["reference"] = reformatted["reference"].str.replace(
        "\\]", "", regex=False
    )

    # Remove duplicate references
    ref_list = []
    for x in reformatted["reference"]:
        y = list(set(x.split("; ")))
        ref_list.append("; ".join(y))
    reformatted["reference"] = ref_list

    fname = os.path.join(outdir, f"{pokay_name}.pokay_database_hits.concise.csv")
    reformatted.drop_duplicates().to_csv(fname, index=None)

    # Print stats
    genes = len(reformatted["amino_acid_consequence"].unique())
    hits = len(reformatted["category"])
    print(f"{genes} new amino acid changes had hits to {hits} categories")
    print(f"See: {fname}")

    return reformatted1
