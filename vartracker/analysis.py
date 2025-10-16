"""
Analysis and plotting functionality for vartracker.

Contains functions for generating plots and analyzing mutation patterns.
"""

import html
import os
import re
import string
from typing import Dict, List, Optional, Sequence, Tuple, Union, cast

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from .constants import REF_GENE_LENGTHS, NSP_LENGTHS, NSPS
from .core import get_logo

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

    matrix_df = pd.DataFrame(records)
    matrix_df = matrix_df.sort_values(["gene_order", "start", "base_label"])
    matrix_df = matrix_df.drop(columns=["gene_order", "start"], errors="ignore")
    matrix_df = matrix_df.drop_duplicates(subset=["base_label"])
    base_label_map = dict(zip(matrix_df["label"], matrix_df["base_label"]))
    canonical_label_map = {
        label: _canonical_label(base_label)
        for label, base_label in base_label_map.items()
    }
    matrix_df = matrix_df.set_index("label")
    matrix_df = matrix_df.drop(columns=["base_label"], errors="ignore")

    matrix_df = matrix_df.loc[matrix_df.max(axis=1) > 0]
    matrix_df = matrix_df.reindex(columns=ordered_samples).fillna(0.0)
    matrix_df.attrs["base_labels"] = base_label_map
    matrix_df.attrs["canonical_labels"] = canonical_label_map

    return matrix_df


def _normalise_base_label(gene: str | None, amino_change: str | None) -> str:
    gene = str(gene or "").strip()
    amino_change = str(amino_change or "").strip()
    if amino_change and amino_change.startswith(f"{gene}:"):
        return amino_change
    if gene and amino_change:
        return f"{gene}:{amino_change}"
    return amino_change or gene


def _slugify_anchor(value: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z]+", "-", value).strip("-").lower()
    if not cleaned:
        cleaned = "entry"
    return f"pokay-{cleaned}"


def _canonical_label(value: str) -> str:
    if not value:
        return ""
    gene, _, aa = str(value).partition(":")
    gene_root = gene.split("_", 1)[0]
    gene_norm = re.sub(r"[^0-9a-z]+", "", gene_root.lower())
    aa_norm = re.sub(r"[^0-9a-z]+", "", aa.lower()) if aa else ""
    return f"{gene_norm}:{aa_norm}" if aa_norm else gene_norm


def _linkify_reference_cell(text: str) -> str:
    if not text:
        return ""

    pattern = re.compile(r"(https?://[^\s;]+)|(10\.\d{4,9}/[^\s;]+)", re.IGNORECASE)
    result: List[str] = []
    last = 0

    for match in pattern.finditer(text):
        result.append(html.escape(text[last : match.start()]))
        token = match.group(0)
        if token.lower().startswith("http"):
            href = token if token.lower().startswith("http") else f"https://{token}"
        else:
            href = f"https://doi.org/{token}"
        link = (
            f'<a href="{html.escape(href)}" target="_blank" rel="noopener">'
            f"{html.escape(token)}</a>"
        )
        result.append(link)
        last = match.end()

    result.append(html.escape(text[last:]))
    return "".join(result)


def _build_pokay_table_html(
    df: Optional[pd.DataFrame],
    anchor_lookup: Dict[str, str],
    table_path: Optional[str],
) -> str:
    if df is None or df.empty:
        return '<p class="pokay-empty">No pokay matches were found for the new mutations.</p>'

    columns = [html.escape(str(col)) for col in df.columns]
    header_cells = "".join(f"<th>{col}</th>" for col in columns)
    header_html = f"<thead><tr>{header_cells}</tr></thead>"

    body_rows: List[str] = []
    primary_rows: set[str] = set()
    for idx, row in df.iterrows():
        base_label = _normalise_base_label(
            row.get("gene"), row.get("amino_acid_consequence")
        )
        canonical_key = _canonical_label(base_label)
        anchor_id = anchor_lookup.get(
            canonical_key, _slugify_anchor(canonical_key or f"{base_label}-{idx}")
        )
        body_cells = []
        for col in df.columns:
            value = row.get(col, "")
            if pd.isna(value):
                value = ""
            if str(col).lower() == "reference":
                rendered = _linkify_reference_cell(str(value))
            else:
                rendered = html.escape(str(value))
            body_cells.append(f"<td>{rendered}</td>")
        if canonical_key not in primary_rows and anchor_id:
            row_attrs = f'id="{anchor_id}" data-anchor="{canonical_key}"'
            primary_rows.add(canonical_key)
        else:
            row_attrs = f'data-anchor="{canonical_key}"'
        body_rows.append(
            f"<tr {row_attrs} class=\"pokay-row\">{''.join(body_cells)}</tr>"
        )

    caption = ""
    if table_path:
        caption = f'<caption class="info-note">Source: {html.escape(os.path.basename(table_path))}</caption>'

    body_html = "<tbody>" + "".join(body_rows) + "</tbody>"
    return f'<table class="pokay-table">{caption}{header_html}{body_html}</table>'


def _write_interactive_heatmap_html(
    matrix: pd.DataFrame,
    sample_names: Sequence[str],
    outdir: str,
    project_name: str,
    pokay_df: Optional[pd.DataFrame],
    pokay_table_path: Optional[str],
    cli_command: Optional[str],
) -> None:
    if matrix.empty:
        return

    x_labels = [str(name) for name in sample_names]
    y_labels = list(matrix.index)
    if not x_labels or not y_labels:
        return

    label_map: Dict[str, str] = matrix.attrs.get("base_labels", {})
    canonical_label_map: Dict[str, str] = matrix.attrs.get("canonical_labels", {})

    if pokay_df is None and pokay_table_path and os.path.exists(pokay_table_path):
        try:
            pokay_df = pd.read_csv(pokay_table_path)
        except Exception:
            pokay_df = None

    highlight_keys: Dict[str, str] = {}
    if pokay_df is not None and not pokay_df.empty:
        for _, prow in pokay_df.iterrows():
            key = _normalise_base_label(
                prow.get("gene"), prow.get("amino_acid_consequence")
            )
            canonical = _canonical_label(key)
            if canonical:
                highlight_keys.setdefault(canonical, _slugify_anchor(canonical))

    cmap = mpl.colormaps.get_cmap("viridis")

    def _frequency_to_color(value: float) -> tuple[str, str]:
        if value is None or np.isnan(value) or value <= 0:
            return "#d9d9d9", ""
        clipped = float(np.clip(value, 0.0, 1.0))
        rgba = cmap(clipped)
        return mpl.colors.to_hex(rgba, keep_alpha=False), f"{clipped:.2f}"

    grid_cells: List[str] = []
    grid_cells.append('<div class="grid-header grid-corner">Variant</div>')
    for sample in x_labels:
        grid_cells.append(
            f'<div class="grid-header grid-sample">{html.escape(sample)}</div>'
        )

    for label in y_labels:
        base_label = label_map.get(label, label.replace("\n", " "))
        canonical_key = canonical_label_map.get(label, _canonical_label(base_label))
        anchor_id = highlight_keys.get(canonical_key)
        safe_label = html.escape(label).replace("\n", "<br>")
        if anchor_id:
            header_html = (
                '<div class="grid-header grid-variant">'
                f'<a href="#pokay-results" class="heatmap-anchor" data-anchor="{canonical_key}">{safe_label}</a>'
                "</div>"
            )
        else:
            header_html = (
                '<div class="grid-header grid-variant">'
                f'<span class="heatmap-label">{safe_label}</span>'
                "</div>"
            )
        grid_cells.append(header_html)

        row_values = matrix.loc[label]
        for sample, value in zip(x_labels, row_values):
            freq = float(value) if value is not None else 0.0
            color, text_value = _frequency_to_color(freq)
            tooltip = html.escape(
                f"Variant: {label.replace(chr(10), ' ')} • "
                f"Sample: {sample} • "
                f"Allele frequency: {freq:.3f}"
            )
            classes = ["cell"]
            if not text_value:
                classes.append("cell-empty")
            cell_html = (
                f'<div class="{" ".join(classes)}" '
                f'style="background-color:{color};" '
                f'data-tooltip="{tooltip}">'
                f'<span class="cell-value">{html.escape(text_value)}</span>'
                "</div>"
            )
            grid_cells.append(cell_html)

    grid_style = f"grid-template-columns: minmax(280px, 320px) repeat({len(x_labels)}, minmax(120px, 1fr));"
    heatmap_grid_html = (
        f'<div class="heatmap-grid" style="{grid_style}">{"".join(grid_cells)}</div>'
    )
    heatmap_scroll_html = f'<div class="heatmap-scroll">{heatmap_grid_html}</div>'

    heatmap_title = (
        f"{project_name}: variant allele frequencies"
        if project_name
        else "Variant allele frequencies"
    )

    pokay_table_html = _build_pokay_table_html(
        pokay_df, highlight_keys, pokay_table_path
    )
    logo_html = f'<pre class="logo">{html.escape(get_logo())}</pre>'
    summary_html = ""
    if cli_command:
        summary_html = (
            '<section class="card summary">'
            "<h1>Workflow summary</h1>"
            f'<p class="summary-text"><code>{html.escape(cli_command)}</code></p>'
            "</section>"
        )

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{html.escape(heatmap_title)}</title>
  <style>
    body {{
      margin: 0;
      font-family: "Inter", "Segoe UI", sans-serif;
      background: #f9fafb;
      color: #111827;
    }}
    main {{
      max-width: 1100px;
      margin: 0 auto;
      padding: 32px 24px 48px;
    }}
    .logo {{
      background: #0f172a;
      color: #e2e8f0;
      padding: 16px;
      border-radius: 12px;
      display: flex;
      justify-content: center;
      align-items: center;
      text-align: center;
      overflow-x: auto;
      font-size: 12px;
      margin-bottom: 24px;
      white-space: pre;
    }}
    .card {{
      background: #ffffff;
      border-radius: 16px;
      box-shadow: 0 12px 24px rgba(15, 23, 42, 0.1);
      padding: 24px;
      margin-bottom: 32px;
    }}
    .card h1 {{
      margin-top: 0;
      font-size: 1.5rem;
      color: #0f172a;
    }}
    .heatmap-scroll {{
      width: 100%;
      overflow-x: auto;
      padding-bottom: 8px;
    }}
    .heatmap-grid {{
      display: grid;
      gap: 6px;
      align-items: stretch;
      width: max(100%, 960px);
    }}
    .grid-header {{
      background: #111827;
      color: #f8fafc;
      padding: 12px 14px;
      border-radius: 12px;
      font-weight: 600;
      display: flex;
      align-items: center;
      justify-content: center;
      text-align: center;
    }}
    .grid-header.grid-variant {{
      justify-content: center;
      text-align: center;
      flex-direction: column;
    }}
    .grid-header.grid-sample {{
      background: #1f2937;
    }}
    .cell {{
      position: relative;
      min-height: 60px;
      border-radius: 12px;
      display: flex;
      align-items: center;
      justify-content: center;
      font-weight: 600;
      color: #0f172a;
      overflow: hidden;
      transition: transform 0.15s ease;
    }}
    .cell-empty {{
      color: #1f2937;
    }}
    .cell:hover {{
      transform: scale(1.03);
      z-index: 2;
    }}
    .cell::after {{
      content: attr(data-tooltip);
      position: absolute;
      left: 50%;
      bottom: 100%;
      transform: translate(-50%, -8px);
      background: rgba(15, 23, 42, 0.92);
      color: #f8fafc;
      padding: 6px 10px;
      border-radius: 6px;
      white-space: nowrap;
      opacity: 0;
      pointer-events: none;
      font-size: 0.75rem;
      transition: opacity 0.15s ease;
    }}
    .cell:hover::after {{
      opacity: 1;
    }}
    .cell-value {{
      pointer-events: none;
      opacity: 0;
      color: #f8fafc;
      transition: opacity 0.15s ease;
    }}
    .cell:hover .cell-value {{
      opacity: 1;
    }}
    .cell-empty .cell-value {{
      display: none;
    }}
    .heatmap-anchor {{
      color: #c0392b;
      text-decoration: none;
      font-weight: 600;
    }}
    .heatmap-anchor:hover {{
      text-decoration: underline;
    }}
    .pokay-table {{
      width: 100%;
      border-collapse: collapse;
      margin-top: 16px;
      font-size: 0.95rem;
    }}
    .pokay-table caption {{
      caption-side: bottom;
      text-align: left;
      padding-top: 8px;
    }}
    .pokay-table th,
    .pokay-table td {{
      border: 1px solid #e5e7eb;
      padding: 8px 12px;
      text-align: left;
      vertical-align: top;
    }}
    .pokay-table th {{
      background: #f3f4f6;
      font-weight: 600;
    }}
    .pokay-table tr:nth-child(even) {{
      background: #f9fafb;
    }}
    .pokay-row.is-active {{
      background: #fff6c7 !important;
      box-shadow: inset 0 0 0 2px #facc15 !important;
    }}
    .pokay-empty {{
      color: #6b7280;
      margin: 0;
    }}
    .info-note {{
      color: #6b7280;
      font-size: 0.9rem;
      margin-top: 12px;
    }}
    .summary-text {{
      margin: 0;
    }}
    .summary-text code {{
      display: block;
      background: #111827;
      color: #f8fafc;
      padding: 12px 14px;
      border-radius: 10px;
      font-family: "Fira Code", "SFMono-Regular", monospace;
      font-size: 0.9rem;
      white-space: pre-wrap;
      word-break: break-word;
      box-shadow: inset 0 0 0 1px rgba(248, 250, 252, 0.1);
    }}
    .table-scroll {{
      width: 100%;
      overflow-x: auto;
      padding-bottom: 8px;
    }}
  </style>
</head>
<body>
  <main>
    {logo_html}
    {summary_html}
    <section class="card">
      <h1>Interactive variant heatmap</h1>
      <p class="info-note">Hover to view exact allele frequencies. Click highlighted variants to jump to supporting pokay annotations.</p>
      {heatmap_scroll_html}
    </section>
    <section class="card" id="pokay-results">
      <h1>Pokay results</h1>
      <div class="table-scroll">{pokay_table_html}</div>
    </section>
  </main>
  <script>
    (function() {{
      const anchors = Array.from(document.querySelectorAll('.heatmap-anchor'));
      const rows = Array.from(document.querySelectorAll('.pokay-row'));
      const clearActive = () => rows.forEach(row => {{
        row.classList.remove('is-active');
        row.style.backgroundColor = '';
        row.style.boxShadow = '';
      }});
      anchors.forEach(anchor => {{
        anchor.addEventListener('click', event => {{
          const target = anchor.getAttribute('data-anchor');
          if (!target) {{ return; }}
          const matches = rows.filter(row => row.getAttribute('data-anchor') === target);
          if (!matches.length) {{ return; }}
          event.preventDefault();
          clearActive();
          matches.forEach(row => {{
            row.classList.add('is-active');
            row.style.backgroundColor = '#fff6c7';
            row.style.boxShadow = 'inset 0 0 0 2px #facc15';
          }});
          matches[0].scrollIntoView({{ behavior: 'smooth', block: 'center' }});
          history.replaceState(null, '', '#pokay-results');
        }});
      }});
    }})();
  </script>
</body>
</html>
"""

    output_path = os.path.join(outdir, "variant_allele_frequency_heatmap.html")
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(html_doc)


def generate_variant_heatmap(
    table: pd.DataFrame,
    sample_names: Sequence[str],
    sample_numbers: Union[Sequence[int], Sequence[str]],
    outdir: str,
    project_name: str,
    min_snv_freq: float,
    min_indel_freq: float,
    gene_lengths: Dict[str, int] | None = None,
    pokay_hits: Optional[pd.DataFrame] = None,
    pokay_table_path: Optional[str] = None,
    cli_command: Optional[str] = None,
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

        try:
            _write_interactive_heatmap_html(
                heatmap_data,
                sample_names,
                outdir,
                project_name,
                pokay_hits,
                pokay_table_path,
                cli_command,
            )
        except Exception as html_exc:  # pragma: no cover - best-effort UX
            print(f"Warning: failed to generate interactive heatmap: {html_exc}")

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
