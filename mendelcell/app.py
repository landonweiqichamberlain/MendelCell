import io
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from mendelcell import list_tissues, run_mendelcell
from mendelcell.report import create_pdf_report, safe_filename


# -----------------------------
# Page setup
# -----------------------------

st.set_page_config(
    page_title="MendelCell",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 MendelCell")
st.subheader("Candidate gene prioritization by tissue-specific single-cell expression")

st.write(
    "Upload Human Protein Atlas single-cell expression files and a candidate gene list. "
    "MendelCell will identify tissue-specific cell types, find candidate genes expressed "
    "above your threshold, and generate tables, plots, and a PDF report."
)


# -----------------------------
# Helper functions
# -----------------------------

@st.cache_data
def read_uploaded_tsv(file_name, file_bytes):
    """
    Read uploaded TSV, TXT, ZIP, or GZ file.
    """
    import tempfile
    import zipfile
    import gzip
    from pathlib import Path

    suffix = Path(file_name).suffix.lower()

    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(file_bytes)
        tmp_path = Path(tmp.name)

    try:
        file_name_lower = file_name.lower()

        if file_name_lower.endswith(".zip"):
            with zipfile.ZipFile(tmp_path) as z:
                file_list = [
                    f for f in z.namelist()
                    if not f.startswith("__MACOSX")
                    and not f.endswith("/")
                    and not Path(f).name.startswith(".")
                ]

                if len(file_list) == 0:
                    raise ValueError("ZIP file is empty.")

                # Use the first TSV/TXT file inside the ZIP
                tsv_files = [
                    f for f in file_list
                    if f.lower().endswith((".tsv", ".txt"))
                ]

                if len(tsv_files) == 0:
                    raise ValueError(
                        f"No TSV/TXT file found inside ZIP. Files found: {file_list}"
                    )

                with z.open(tsv_files[0]) as f:
                    return pd.read_csv(f, sep="\t", low_memory=False)

        if file_name_lower.endswith(".gz"):
            with gzip.open(tmp_path, "rt") as f:
                return pd.read_csv(f, sep="\t", low_memory=False)

        return pd.read_csv(tmp_path, sep="\t", low_memory=False)

    except Exception as e:
        raise ValueError(f"Could not read {file_name}: {e}")

    finally:
        tmp_path.unlink(missing_ok=True)


def make_gene_count_plot(results):
    fig, ax = plt.subplots(figsize=(11, 6))

    plot_df = results.plot_df

    ax.bar(plot_df["Cell type"], plot_df["Gene count"])
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Number of candidate genes")
    ax.set_title(
        f"Candidate genes expressed in {results.selected_tissue}-specific cell types"
    )

    ax.tick_params(axis="x", rotation=45)

    for label in ax.get_xticklabels():
        label.set_ha("right")

    fig.tight_layout()
    return fig


def make_ncpm_plot(results, cell_type):
    cell_df = results.ncpm_df[results.ncpm_df["Cell type"] == cell_type]
    cell_df = cell_df.sort_values("nCPM", ascending=False)

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.bar(cell_df["Gene name"], cell_df["nCPM"])
    ax.set_xlabel("Gene")
    ax.set_ylabel("Mean nCPM")
    ax.set_title(f"Mean nCPM of candidate genes in {cell_type}")

    ax.tick_params(axis="x", rotation=45)

    for label in ax.get_xticklabels():
        label.set_ha("right")

    fig.tight_layout()
    return fig


# -----------------------------
# Sidebar upload widgets
# -----------------------------

st.sidebar.header("1. Upload input files")

cluster_file = st.sidebar.file_uploader(
    "Upload rna_single_cell_cluster.tsv or .zip",
    type=["tsv", "txt", "zip"]
)

hpa_file = st.sidebar.file_uploader(
    "Upload rna_single_cell_type.tsv or .zip",
    type=["tsv", "txt", "zip"]
)

gene_file = st.sidebar.file_uploader(
    "Upload candidate gene list TSV",
    type=["tsv", "txt"]
)

if cluster_file is None or hpa_file is None or gene_file is None:
    st.info("Upload all three files to begin.")

    st.markdown(
        """
        Required files:

        1. `rna_single_cell_cluster.tsv` or `.zip`
        2. `rna_single_cell_type.tsv` or `.zip`
        3. Candidate gene list TSV with a column named `Gene Symbol`
        """
    )

    st.stop()


# -----------------------------
# Read uploaded files
# -----------------------------

try:
    st.write("Reading cluster file:", cluster_file.name)
    clusters = read_uploaded_tsv(cluster_file.name, cluster_file.getvalue())
    st.success(f"Cluster file loaded: {clusters.shape[0]} rows, {clusters.shape[1]} columns")

    st.write("Reading HPA cell-type file:", hpa_file.name)
    hpa = read_uploaded_tsv(hpa_file.name, hpa_file.getvalue())
    st.success(f"HPA file loaded: {hpa.shape[0]} rows, {hpa.shape[1]} columns")

    st.write("Reading gene list file:", gene_file.name)
    gene_table = read_uploaded_tsv(gene_file.name, gene_file.getvalue())
    st.success(f"Gene list loaded: {gene_table.shape[0]} rows, {gene_table.shape[1]} columns")

except Exception as e:
    st.error(f"Could not read uploaded files: {e}")
    st.exception(e)
    st.stop()




# -----------------------------
# Sidebar analysis settings
# -----------------------------

try:
    valid_tissues = list_tissues(clusters)

except Exception as e:
    st.error(f"Could not find tissue names in cluster file: {e}")
    st.stop()

st.sidebar.header("2. Choose settings")

selected_tissue = st.sidebar.selectbox(
    "Select tissue",
    valid_tissues
)

threshold = st.sidebar.number_input(
    "Expression threshold",
    min_value=0.0,
    value=1.0,
    step=0.5
)

run_button = st.sidebar.button("Run MendelCell analysis")


# -----------------------------
# Run analysis
# -----------------------------

if not run_button:
    st.info("Choose a tissue and threshold, then click **Run MendelCell analysis**.")
    st.stop()

try:
    results = run_mendelcell(
        clusters=clusters,
        hpa=hpa,
        gene_table=gene_table,
        tissue=selected_tissue,
        threshold=threshold
    )

except Exception as e:
    st.error(f"MendelCell analysis failed: {e}")
    st.stop()


# -----------------------------
# Display results
# -----------------------------

st.success("Analysis complete.")

col1, col2, col3, col4 = st.columns(4)

col1.metric("Input genes", len(results.gene_symbols))
col2.metric("Tissue-specific cell types", len(results.unique_cells))
col3.metric("Genes passing threshold", results.filtered["Gene name"].nunique())
col4.metric(
    "Gene-cell pairs",
    len(results.filtered[["Cell type", "Gene name"]].drop_duplicates())
)

st.header("Cell types unique to selected tissue")
st.dataframe(results.unique_to_tissue, use_container_width=True)

st.header("Candidate gene count per cell type")
st.dataframe(results.cell_count_df, use_container_width=True)

if not results.plot_df.empty:
    fig = make_gene_count_plot(results)
    st.pyplot(fig)
    plt.close(fig)

st.header("Candidate genes found in each cell type")
st.dataframe(results.genes_in_cell_df, use_container_width=True)

st.header("Filtered candidate genes")
st.dataframe(results.filtered_report, use_container_width=True)

st.header("Mean nCPM values")
st.dataframe(results.ncpm_df, use_container_width=True)

st.header("nCPM plots by cell type")

for cell_type in results.ncpm_df["Cell type"].unique():
    fig = make_ncpm_plot(results, cell_type)
    st.pyplot(fig)
    plt.close(fig)


# -----------------------------
# Create downloadable outputs
# -----------------------------

st.header("Download outputs")

unique_tsv = results.unique_to_tissue.to_csv(sep="\t", index=False)
filtered_tsv = results.filtered.to_csv(sep="\t", index=False)
ncpm_tsv = results.ncpm_df.to_csv(sep="\t", index=False)

safe_tissue = safe_filename(results.selected_tissue)
pdf_filename = f"MendelCell_report_{safe_tissue}.pdf"

with tempfile.TemporaryDirectory() as tmpdir:
    pdf_path = Path(tmpdir) / pdf_filename
    create_pdf_report(results, pdf_path)

    with open(pdf_path, "rb") as f:
        pdf_bytes = f.read()

st.download_button(
    label="Download PDF report",
    data=pdf_bytes,
    file_name=pdf_filename,
    mime="application/pdf"
)

st.download_button(
    label="Download unique cell types TSV",
    data=unique_tsv,
    file_name="unique_cell_types.tsv",
    mime="text/tab-separated-values"
)

st.download_button(
    label="Download filtered candidate genes TSV",
    data=filtered_tsv,
    file_name="filtered_candidate_genes.tsv",
    mime="text/tab-separated-values"
)

st.download_button(
    label="Download nCPM table TSV",
    data=ncpm_tsv,
    file_name="candidate_gene_ncpm_by_cell_type.tsv",
    mime="text/tab-separated-values"
)