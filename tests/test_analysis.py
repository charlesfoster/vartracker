"""Tests for analysis helpers."""

from __future__ import annotations

import os
import string

import matplotlib.pyplot as plt
import pandas as pd
import pytest
import seaborn as sns

from vartracker.analysis import (
    _heatmap_figure_size,
    search_literature,
    _prepare_variant_heatmap_matrix,
    process_joint_variants,
    generate_variant_heatmap,
    generate_gene_table,
    plot_gene_table,
    parse_plot_genes_arg,
    select_genes_for_plot,
)


def _gene_table_rows(gene_table, gene):
    rows = gene_table[gene_table["gene"] == gene]
    return dict(zip(rows["type"], rows["number"]))


def test_generate_gene_table_keeps_single_contig_genes_unqualified():
    """Baseline: single-contig data (e.g. viral references) is unaffected."""
    table = pd.DataFrame(
        {
            "gene": ["S", "S", "N"],
            "chrom": ["NC_045512.2", "NC_045512.2", "NC_045512.2"],
            "type_of_change": ["missense", "missense", "synonymous"],
            "presence_absence": ["NY", "NY", "YY"],
        }
    )

    gene_table = generate_gene_table(table, gene_lengths={"S": 3822, "N": 1260})

    assert set(gene_table["gene"]) == {"S", "N"}
    assert _gene_table_rows(gene_table, "S")["total"] == 2


def test_generate_gene_table_splits_colliding_gene_names_across_contigs():
    """Same gene name on chromosome + plasmid must not be summed together."""
    table = pd.DataFrame(
        {
            "gene": ["repA", "repA", "dnaA"],
            "chrom": ["chrom1", "plasmid1", "chrom1"],
            "type_of_change": ["missense", "missense", "missense"],
            "presence_absence": ["NY", "NY", "NY"],
        }
    )

    gene_table = generate_gene_table(
        table,
        gene_lengths={
            "dnaA": 300,
            "repA (chrom1)": 201,
            "repA (plasmid1)": 150,
        },
    )

    genes_present = set(gene_table["gene"])
    assert "repA (chrom1)" in genes_present
    assert "repA (plasmid1)" in genes_present
    assert "repA" not in genes_present

    # Each contig's copy of repA should retain its own count, not a merged total.
    assert _gene_table_rows(gene_table, "repA (chrom1)")["total"] == 1
    assert _gene_table_rows(gene_table, "repA (plasmid1)")["total"] == 1


def test_generate_gene_table_ambiguous_genes_param_prevents_data_loss():
    """A gene flagged ambiguous by the annotation, but with variants on only
    one contig in this particular table, must still surface under its
    qualified label so it matches the gene-length scaffold and isn't dropped
    by the scaffold merge."""
    table = pd.DataFrame(
        {
            "gene": ["repA"],
            "chrom": ["chrom1"],
            "type_of_change": ["missense"],
            "presence_absence": ["NY"],
        }
    )

    gene_table = generate_gene_table(
        table,
        gene_lengths={
            "repA (chrom1)": 201,
            "repA (plasmid1)": 150,
        },
        ambiguous_genes={"repA"},
    )

    genes_present = set(gene_table["gene"])
    assert "repA (chrom1)" in genes_present
    assert "repA" not in genes_present
    assert _gene_table_rows(gene_table, "repA (chrom1)")["total"] == 1
    # The scaffold still surfaces the other contig's copy with a zero count.
    assert _gene_table_rows(gene_table, "repA (plasmid1)")["total"] == 0


def _make_variant_row(gene, type_of_change, presence_absence, chrom="chrom1"):
    return {
        "gene": gene,
        "chrom": chrom,
        "type_of_change": type_of_change,
        "presence_absence": presence_absence,
    }


def test_select_genes_for_plot_keeps_everything_when_annotation_under_cap():
    """When the whole annotated gene set already fits under the cap (e.g. a
    viral reference), nothing is filtered - this is what keeps SARS-CoV-2
    figures byte-identical regardless of the new default cap."""
    table = pd.DataFrame(
        [
            _make_variant_row("S", "missense", "NY"),
            _make_variant_row("N", "synonymous", "YY"),
        ]
    )
    gene_table = generate_gene_table(
        table, gene_lengths={"S": 3822, "N": 1260, "E": 228}
    )

    plot_table, subtitle = select_genes_for_plot(gene_table, max_plot_genes=30)

    assert subtitle is None
    pd.testing.assert_frame_equal(
        plot_table.reset_index(drop=True), gene_table.reset_index(drop=True)
    )


def test_select_genes_for_plot_ranks_by_new_mutations_with_total_tiebreak():
    """Ranking must use newly emerged variants first, falling back to total
    variants only to break ties - not the other way around."""
    table = pd.DataFrame(
        [
            # Gene A: 3 new mutations, 5 total.
            _make_variant_row("A", "missense", "NY"),
            _make_variant_row("A", "missense", "NY"),
            _make_variant_row("A", "missense", "NY"),
            _make_variant_row("A", "synonymous", "YY"),
            _make_variant_row("A", "synonymous", "YY"),
            # Gene B: 3 new mutations (tied with A), only 2 total.
            _make_variant_row("B", "missense", "NY"),
            _make_variant_row("B", "missense", "NY"),
            _make_variant_row("B", "missense", "NY"),
            # Gene C: fewer new mutations than A/B, but more total variants.
            _make_variant_row("C", "missense", "NY"),
            *[_make_variant_row("C", "synonymous", "YY") for _ in range(9)],
            # Gene D: no new mutations at all.
            _make_variant_row("D", "synonymous", "YY"),
        ]
    )
    gene_table = generate_gene_table(
        table, gene_lengths={"A": 100, "B": 100, "C": 100, "D": 100, "E": 100}
    )

    plot_table, subtitle = select_genes_for_plot(gene_table, max_plot_genes=2)

    assert subtitle == "top 2 of 4 genes with variants"
    assert set(plot_table["gene"]) == {"A", "B"}


def test_select_genes_for_plot_omits_subtitle_when_variant_genes_fit_cap():
    """Even if the annotated genome is large, if the genes that actually
    carry variants already fit under the cap, nothing was truncated and the
    subtitle must be omitted (e.g. never print "top 30 of 12")."""
    table = pd.DataFrame(
        [
            _make_variant_row("A", "missense", "NY"),
            _make_variant_row("B", "missense", "NY"),
        ]
    )
    gene_lengths = {"A": 100, "B": 100}
    gene_lengths.update({f"unused{i}": 100 for i in range(20)})
    gene_table = generate_gene_table(table, gene_lengths=gene_lengths)

    plot_table, subtitle = select_genes_for_plot(gene_table, max_plot_genes=10)

    assert subtitle is None
    assert set(plot_table["gene"]) == {"A", "B"}


def test_select_genes_for_plot_cap_zero_or_negative_raises():
    table = pd.DataFrame([_make_variant_row("A", "missense", "NY")])
    gene_table = generate_gene_table(table, gene_lengths={"A": 100})

    with pytest.raises(ValueError):
        select_genes_for_plot(gene_table, max_plot_genes=0)

    with pytest.raises(ValueError):
        select_genes_for_plot(gene_table, max_plot_genes=-5)


def test_select_genes_for_plot_explicit_list_warns_and_drops_invalid(capsys):
    table = pd.DataFrame(
        [
            _make_variant_row("X", "missense", "NY"),
            _make_variant_row("Y", "missense", "NY"),
        ]
    )
    gene_lengths = {"X": 100, "Y": 100, "Z": 100}  # Z is annotated but has no variants
    gene_table = generate_gene_table(table, gene_lengths=gene_lengths)

    plot_table, subtitle = select_genes_for_plot(
        gene_table, plot_genes=["X", "Z", "FAKE"]
    )

    assert subtitle is None
    assert set(plot_table["gene"]) == {"X"}
    captured = capsys.readouterr()
    assert "Z" in captured.out and "no variants" in captured.out
    assert "FAKE" in captured.out and "reference annotation" in captured.out


def test_select_genes_for_plot_explicit_list_all_invalid_raises():
    table = pd.DataFrame([_make_variant_row("X", "missense", "NY")])
    gene_table = generate_gene_table(table, gene_lengths={"X": 100})

    with pytest.raises(ValueError):
        select_genes_for_plot(gene_table, plot_genes=["FAKE1", "FAKE2"])


def test_select_genes_for_plot_explicit_list_overrides_max_plot_genes():
    table = pd.DataFrame(
        [
            _make_variant_row("A", "missense", "NY"),
            _make_variant_row("B", "missense", "NY"),
            _make_variant_row("C", "missense", "NY"),
        ]
    )
    gene_table = generate_gene_table(table, gene_lengths={"A": 100, "B": 100, "C": 100})

    plot_table, subtitle = select_genes_for_plot(
        gene_table, max_plot_genes=1, plot_genes=["A", "B"]
    )

    assert subtitle is None
    assert set(plot_table["gene"]) == {"A", "B"}


def test_parse_plot_genes_arg_comma_separated_deduplicates_and_preserves_order():
    genes = parse_plot_genes_arg("geneA, geneB,geneA, geneC")
    assert genes == ["geneA", "geneB", "geneC"]


def test_parse_plot_genes_arg_reads_from_file(tmp_path):
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("geneA\n\n# a comment\ngeneB\ngeneA\n")

    genes = parse_plot_genes_arg(str(gene_file))

    assert genes == ["geneA", "geneB"]


def _golden_plot_gene_table(gene_table, pname, outdir):
    """Frozen copy of `plot_gene_table` as it existed before the
    `--max-plot-genes`/`--plot-genes` options were added. Used only as a
    ground truth for the byte-identical regression test below."""
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
            label.set_rotation(90)

    for ax, title in zip(g.fig.axes, list(gene_table["type"].unique())):
        if pd.isna(title):
            title_str = "None"
        else:
            title_str = str(title) if not isinstance(title, str) else title
        ax.set_title(string.capwords(title_str.replace("_", " ")))

    plt.savefig(
        os.path.join(outdir, "mutations_per_gene.pdf"), dpi=300, bbox_inches="tight"
    )
    plt.close()


def test_plot_gene_table_default_matches_pre_cap_baseline(tmp_path, monkeypatch):
    """With the default cap (30) and a SARS-CoV-2-sized (12 gene) reference,
    output must be byte-identical to the pre-existing plotting behaviour, so
    published manuscript figures don't change."""
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")

    gene_names = [f"gene{i}" for i in range(12)]
    rows = []
    for i, gene in enumerate(gene_names[:10]):
        rows.append(_make_variant_row(gene, "missense", "NY" if i % 2 == 0 else "YY"))
    table = pd.DataFrame(rows)
    gene_lengths = {gene: 1000 for gene in gene_names}
    gene_table = generate_gene_table(table, gene_lengths=gene_lengths)

    golden_dir = tmp_path / "golden"
    new_dir = tmp_path / "new"
    golden_dir.mkdir()
    new_dir.mkdir()

    _golden_plot_gene_table(gene_table.copy(deep=True), "Baseline", str(golden_dir))
    plot_gene_table(gene_table.copy(deep=True), "Baseline", str(new_dir))

    golden_bytes = (golden_dir / "mutations_per_gene.pdf").read_bytes()
    new_bytes = (new_dir / "mutations_per_gene.pdf").read_bytes()
    assert golden_bytes == new_bytes


def test_search_literature_handles_nullable_boolean_masks(tmp_path):
    table = pd.DataFrame(
        {
            "gene": ["S"],
            "amino_acid_consequence": ["S:N501Y"],
            "variant_status": ["new"],
        }
    )

    literature = pd.DataFrame(
        {
            "gene": pd.Series(["S", pd.NA], dtype="string"),
            "mutation": pd.Series(["S:N501Y", ""], dtype="string"),
            "category": ["Functional impact", ""],
            "information": ["Info", ""],
            "reference": ["PMID123", ""],
        }
    )

    result = search_literature(table, literature, tmp_path, "sample")

    assert not result.empty
    assert (tmp_path / "sample.literature_database_hits.full.csv").exists()
    assert (tmp_path / "sample.literature_database_hits.concise.csv").exists()


def test_prepare_variant_heatmap_matrix_orders_variants_by_genome():
    table = pd.DataFrame(
        [
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "809T",
                "nsp_aa_change": "",
                "type_of_change": "synonymous",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.5",
                "samples": "P0 / P1",
                "variant": "C1059T",
                "start": 1059,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 1.0",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {  # duplicate row should be ignored
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 1.0",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {  # should be dropped: below SNV threshold
                "gene": "S",
                "amino_acid_consequence": "A10T",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.1",
                "samples": "P0 / P1",
                "variant": "G21500A",
                "start": 21500,
            },
            {  # should be dropped: indel below threshold
                "gene": "S",
                "amino_acid_consequence": "144del",
                "nsp_aa_change": "",
                "type_of_change": "inframe_deletion",
                "type_of_variant": "indel",
                "alt_freq": "0.05 / 0.05",
                "samples": "P0 / P1",
                "variant": "CT21991C",
                "start": 21991,
            },
            {
                "gene": "5' UTR",
                "amino_acid_consequence": "",
                "nsp_aa_change": "",
                "type_of_change": "None",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.9",
                "samples": "P0 / P1",
                "variant": "C200T",
                "start": 200,
            },
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "",
                "nsp_aa_change": "nsp2:VLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYIVEELKAA",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.8",
                "samples": "P0 / P1",
                "variant": "G1946GT",
                "start": 1946,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(table, ["P0", "P1"], 0.2, 0.3)

    long_aa = "VLQKAAITILDGISQYSLRLIDAMMFTSDLATNNLVVMAYIVEELKAA"
    head = long_aa[:10]
    tail = long_aa[-10:]
    middle = len(long_aa) - 20

    expected_long_label = f"nsp2:{head}+{middle}{tail}\n(G1946GT)"

    expected_index = [
        "nsp2:T629=\n(C1059T)",
        expected_long_label,
        "S:D215G\n(A22206G)",
    ]

    assert list(matrix.index) == expected_index
    assert matrix.loc["nsp2:T629=\n(C1059T)", "P0"] == 0.0
    assert matrix.loc["nsp2:T629=\n(C1059T)", "P1"] == 0.5
    assert matrix.loc[expected_long_label, "P1"] == 0.8
    assert matrix.loc["S:D215G\n(A22206G)", "P1"] == 1.0


def test_prepare_variant_heatmap_matrix_normalises_starred_synonymous_label():
    table = pd.DataFrame(
        [
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "924F",
                "nsp_aa_change": "",
                "type_of_change": "*synonymous",
                "type_of_variant": "snp",
                "alt_freq": "0.5 / 0.0",
                "samples": "P0 / P1",
                "variant": "C3037T",
                "start": 3037,
            },
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "924F",
                "nsp_aa_change": "",
                "type_of_change": "synonymous",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.6",
                "samples": "P0 / P1",
                "variant": "C3037T",
                "start": 3037,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(table, ["P0", "P1"], 0.2, 0.3)

    assert list(matrix.index) == ["nsp3_PLpro:F106=\n(C3037T)"]
    assert matrix.loc["nsp3_PLpro:F106=\n(C3037T)", "P0"] == 0.5
    assert matrix.loc["nsp3_PLpro:F106=\n(C3037T)", "P1"] == 0.6


def test_prepare_variant_heatmap_matrix_repairs_stale_stop_gained_nsp_label():
    table = pd.DataFrame(
        [
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "L889*",
                "nsp_aa_change": "nsp3_PLpro:71L",
                "type_of_change": "stop_gained",
                "type_of_variant": "snp",
                "alt_freq": "0.052",
                "samples": "P0",
                "variant": "T2931A",
                "start": 2931,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(table, ["P0"], 0.0, 0.0)

    assert list(matrix.index) == ["nsp3_PLpro:L71*\n(T2931A)"]


def test_prepare_variant_heatmap_matrix_excludes_selected_consequence_types():
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.7",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "T716T",
                "nsp_aa_change": "",
                "type_of_change": "synonymous",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.8",
                "samples": "P0 / P1",
                "variant": "C23403T",
                "start": 23403,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "145del",
                "nsp_aa_change": "",
                "type_of_change": "frameshift",
                "type_of_variant": "indel",
                "alt_freq": "0.0 / 0.9",
                "samples": "P0 / P1",
                "variant": "A22029-",
                "start": 22029,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1"],
        0.2,
        0.2,
        excluded_consequence_types=["synonymous", "frameshift"],
    )

    assert list(matrix.index) == ["S:D215G\n(A22206G)"]


def test_prepare_variant_heatmap_matrix_excludes_wildcard_consequence_types():
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "joint_frameshift",
                "type_of_variant": "indel",
                "alt_freq": "0.0 / 0.7",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "N501Y",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.8",
                "samples": "P0 / P1",
                "variant": "A23063T",
                "start": 23063,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1"],
        0.2,
        0.2,
        excluded_consequence_types=["*frameshift*"],
    )

    assert list(matrix.index) == ["S:N501Y\n(A23063T)"]


def test_prepare_variant_heatmap_matrix_applies_extended_filters():
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "all_samples_pass_qc": True,
                "proportion_samples_passing_qc": 1.0,
                "presence_absence": "N / Y / Y",
                "alt_freq": "0.0 / 0.6 / 0.7",
                "variant_site_depth": "0 / 120 / 125",
                "samples": "P0 / P1 / P2",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "145fs",
                "nsp_aa_change": "",
                "type_of_change": "joint_frameshift",
                "type_of_variant": "indel",
                "variant_status": "new",
                "persistence_status": "new_persistent",
                "all_samples_pass_qc": True,
                "proportion_samples_passing_qc": 0.33,
                "presence_absence": "N / N / Y",
                "alt_freq": "0.0 / 0.0 / 0.8",
                "variant_site_depth": "0 / 0 / 130",
                "samples": "P0 / P1 / P2",
                "variant": "A22029-",
                "start": 22029,
            },
            {
                "gene": "N",
                "amino_acid_consequence": "R203K",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "variant_status": "original",
                "persistence_status": "original_retained",
                "all_samples_pass_qc": False,
                "proportion_samples_passing_qc": 0.67,
                "presence_absence": "Y / Y / Y",
                "alt_freq": "0.5 / 0.5 / 0.5",
                "variant_site_depth": "110 / 115 / 120",
                "samples": "P0 / P1 / P2",
                "variant": "G28881A",
                "start": 28881,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1", "P2"],
        0.2,
        0.2,
        only_new=True,
        gene_include=["S"],
        variant_type_include=["snp"],
        qc_include=["true"],
        min_prop_passing_qc=0.9,
        min_persistence=2,
        min_max_af=0.6,
        sample_subset=["P1", "P2"],
        hide_singletons=True,
        min_depth=100,
    )

    assert list(matrix.index) == ["S:D215G\n(A22206G)"]


def test_prepare_variant_heatmap_matrix_only_persistent_includes_new_intermittent():
    """--only-persistent must keep new_intermittent alongside new_persistent,
    since both reached the final timepoint - only the path there differs."""
    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "persistence_status": "new_persistent",
                "alt_freq": "0.0 / 0.6 / 0.7",
                "samples": "P0 / P1 / P2",
                "variant": "A22206G",
                "start": 22206,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "E484K",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "persistence_status": "new_intermittent",
                "alt_freq": "0.0 / 0.6 / 0.7",
                "samples": "P0 / P1 / P2",
                "variant": "G23012A",
                "start": 23012,
            },
            {
                "gene": "S",
                "amino_acid_consequence": "N501Y",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "persistence_status": "new_transient",
                "alt_freq": "0.0 / 0.6 / 0.0",
                "samples": "P0 / P1 / P2",
                "variant": "A23063T",
                "start": 23063,
            },
        ]
    )

    matrix = _prepare_variant_heatmap_matrix(
        table,
        ["P0", "P1", "P2"],
        0.2,
        0.2,
        only_persistent=True,
    )

    assert set(matrix.index) == {"S:D215G\n(A22206G)", "S:E484K\n(G23012A)"}


def test_process_joint_variants_only_adds_single_joint_prefix(tmp_path):
    csv_path = tmp_path / "results.csv"
    pd.DataFrame(
        [
            {
                "start": 100,
                "gene": "S",
                "amino_acid_consequence": "N501Y",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "c.1A>T",
                "bcsq_aa_notation": "p.N501Y",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "joint_*frameshift",
            },
            {
                "start": 101,
                "gene": "",
                "amino_acid_consequence": "",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "",
                "bcsq_aa_notation": "@100",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "frameshift",
            },
        ]
    ).to_csv(csv_path, index=False)

    result = process_joint_variants(str(csv_path))

    assert result.loc[0, "type_of_change"] == "joint_frameshift"
    assert result.loc[1, "type_of_change"] == "joint_frameshift"


def test_process_joint_variants_matches_main_row_by_presence_pattern(tmp_path):
    csv_path = tmp_path / "results.csv"
    pd.DataFrame(
        [
            {
                "start": 100,
                "gene": "S",
                "amino_acid_consequence": "K2A",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "4A>G+5A>C",
                "bcsq_aa_notation": "2K>2A",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "missense",
                "presence_absence": "Y / N",
            },
            {
                "start": 100,
                "gene": "S",
                "amino_acid_consequence": "K2V",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "4A>T+5A>G",
                "bcsq_aa_notation": "2K>2V",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "missense",
                "presence_absence": "N / Y",
            },
            {
                "start": 101,
                "gene": "",
                "amino_acid_consequence": "",
                "nsp_aa_change": "",
                "bcsq_nt_notation": "",
                "bcsq_aa_notation": "@100",
                "aa1_total_properties": "",
                "aa2_total_properties": "",
                "aa1_unique_properties": "",
                "aa2_unique_properties": "",
                "aa1_weight": "",
                "aa2_weight": "",
                "weight_difference": "",
                "type_of_change": "@100",
                "presence_absence": "N / Y",
            },
        ]
    ).to_csv(csv_path, index=False)

    result = process_joint_variants(str(csv_path))

    assert result.loc[2, "amino_acid_consequence"] == "K2V"
    assert result.loc[2, "type_of_change"] == "joint_missense"


def test_process_joint_variants_is_order_invariant_for_overlapping_gene_rows(tmp_path):
    shared_rows = [
        {
            "start": 25470,
            "gene": "ORF3a",
            "amino_acid_consequence": "ORF3a:A10V",
            "nsp_aa_change": "",
            "bcsq_nt_notation": "c.30C>T",
            "bcsq_aa_notation": "p.A10V",
            "aa1_total_properties": "",
            "aa2_total_properties": "",
            "aa1_unique_properties": "",
            "aa2_unique_properties": "",
            "aa1_weight": "",
            "aa2_weight": "",
            "weight_difference": "",
            "type_of_change": "missense",
            "presence_absence": "Y / N / Y",
        },
        {
            "start": 25470,
            "gene": "ORF3c",
            "amino_acid_consequence": "ORF3c:M5I",
            "nsp_aa_change": "",
            "bcsq_nt_notation": "c.15G>A",
            "bcsq_aa_notation": "p.M5I",
            "aa1_total_properties": "",
            "aa2_total_properties": "",
            "aa1_unique_properties": "",
            "aa2_unique_properties": "",
            "aa1_weight": "",
            "aa2_weight": "",
            "weight_difference": "",
            "type_of_change": "missense",
            "presence_absence": "Y / N / Y",
        },
        {
            "start": 25471,
            "gene": "",
            "amino_acid_consequence": "",
            "nsp_aa_change": "",
            "bcsq_nt_notation": "",
            "bcsq_aa_notation": "@25470",
            "aa1_total_properties": "",
            "aa2_total_properties": "",
            "aa1_unique_properties": "",
            "aa2_unique_properties": "",
            "aa1_weight": "",
            "aa2_weight": "",
            "weight_difference": "",
            "type_of_change": "@25470",
            "presence_absence": "Y / N / Y",
        },
    ]

    first_csv = tmp_path / "overlap_first.csv"
    second_csv = tmp_path / "overlap_second.csv"
    pd.DataFrame(shared_rows).to_csv(first_csv, index=False)
    pd.DataFrame([shared_rows[1], shared_rows[0], shared_rows[2]]).to_csv(
        second_csv, index=False
    )

    first_result = process_joint_variants(str(first_csv))
    second_result = process_joint_variants(str(second_csv))

    first_joint = first_result.loc[2, ["gene", "amino_acid_consequence"]].to_dict()
    second_joint = second_result.loc[2, ["gene", "amino_acid_consequence"]].to_dict()

    assert first_joint == second_joint


def test_generate_variant_heatmap_creates_interactive_html(tmp_path, monkeypatch):
    mpl_dir = tmp_path / "mpl"
    mpl_dir.mkdir()
    monkeypatch.setenv("MPLCONFIGDIR", str(mpl_dir))

    table = pd.DataFrame(
        [
            {
                "gene": "S",
                "amino_acid_consequence": "D215G",
                "nsp_aa_change": "",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.0 / 0.5",
                "per_sample_variant_qc": "P / F",
                "samples": "P0 / P1",
                "variant": "A22206G",
                "start": 22206,
                "variant_status": "new",
                "presence_absence": "N / Y",
            },
            {
                "gene": "ORF1ab",
                "amino_acid_consequence": "",
                "nsp_aa_change": "nsp6_TM:L37F",
                "type_of_change": "missense",
                "type_of_variant": "snp",
                "alt_freq": "0.2 / 0.6",
                "per_sample_variant_qc": "P / P",
                "samples": "P0 / P1",
                "variant": "C21575T",
                "start": 21575,
                "variant_status": "new",
                "presence_absence": "N / Y",
            },
        ]
    )

    literature_hits = pd.DataFrame(
        {
            "gene": ["S", "nsp6"],
            "amino_acid_consequence": ["D215G", "L37F"],
            "database_mutation_string": ["S:D215G", "nsp6:L37F"],
            "category": ["Functional impact", "Functional impact"],
            "prior_information": ["Info", "Additional"],
            "reference": ["PMID123; 10.1000/xyz123", "10.2000/abc456"],
        }
    )

    literature_csv = tmp_path / "sample.literature_database_hits.full.csv"
    literature_hits.to_csv(literature_csv, index=False)

    command_str = "vartracker vcf input.csv --outdir results"

    generate_variant_heatmap(
        table,
        ["P0", "P1"],
        [0, 1],
        str(tmp_path),
        "Example",
        0.05,
        0.05,
        literature_hits=literature_hits,
        literature_table_path=str(literature_csv),
        cli_command=command_str,
    )

    pdf_path = tmp_path / "variant_allele_frequency_heatmap.pdf"
    html_path = tmp_path / "variant_allele_frequency_heatmap.html"

    assert pdf_path.exists()
    assert html_path.exists()

    content = html_path.read_text(encoding="utf-8")
    assert "Interactive variant heatmap" in content
    assert "Workflow summary" in content
    assert command_str in content
    assert "heatmap-grid" in content
    assert "heatmap-scroll" in content
    assert "Literature results" in content
    assert "table-scroll" in content
    assert "heatmap-anchor" in content
    assert "cell-qc-fail" in content
    assert "AF=0.50, QC=FAIL" in content
    assert 'data-anchor="s:d215g' in content
    assert 'data-anchor="nsp6:l37f' in content
    assert ">10.1000/xyz123</a>" in content
    assert ".cell:hover .cell-value" in content
    assert "clearActive" in content


def test_heatmap_figure_size_enforces_minimum_row_height():
    width, height = _heatmap_figure_size(6, 2)

    assert width == 4.0
    assert height == 4.0

    width, height = _heatmap_figure_size(8, 2)

    assert width == 4.0
    assert height == pytest.approx(5.2)
