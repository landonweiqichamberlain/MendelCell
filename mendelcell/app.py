from pathlib import Path
import base64
import io
import tempfile

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from mendelcell import list_tissues, run_mendelcell
from mendelcell.report import create_pdf_report, safe_filename


# -----------------------------
# File paths
# -----------------------------

ROOT_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT_DIR / "data"
EXAMPLES_DIR = ROOT_DIR / "examples"

CLUSTER_REFERENCE = DATA_DIR / "mendelcell_clusters_reference.parquet"
CELLTYPE_REFERENCE = DATA_DIR / "mendelcell_celltype_reference.parquet"
EXAMPLE_GENE_LIST = EXAMPLES_DIR / "example_gene_list.tsv"


# -----------------------------
# Page setup
# -----------------------------

st.set_page_config(
    page_title="MendelCell",
    page_icon="🧬",
    layout="wide",
)

st.markdown(
    """
    <style>
    html {
        overflow-y: scroll;
        scrollbar-gutter: stable;
    }

    body {
        overflow-x: hidden;
    }

    [data-testid="stAppViewContainer"] {
        overflow-x: hidden;
    }

    .block-container {
        max-width: 1400px;
        padding-left: 1rem;
        padding-right: 1rem;
        padding-top: 2rem;
    }

    div[data-testid="stDataFrame"] {
        max-width: 1200px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 MendelCell")
st.subheader("Candidate gene prioritization by tissue-specific single-cell expression")

st.write(
    "MendelCell uses a preprocessed Human Protein Atlas single-cell reference. "
    "Upload a candidate gene list, choose a tissue, choose expression settings, "
    "and generate ranked tables, plots, and a PDF report."
)

st.warning(
    "MendelCell is intended for research and exploratory analysis only. "
    "It is not a diagnostic tool and should not be used to make clinical decisions."
)


# -----------------------------
# Helper functions
# -----------------------------

def show_dataframe_with_1_index(df, height=400, width=1200):
    """
    Display a dataframe in Streamlit with row numbering starting at 1.
    """
    display_df = df.copy()
    display_df.index = range(1, len(display_df) + 1)

    st.dataframe(
        display_df,
        width=width,
        height=height,
    )


def show_matplotlib_svg(fig, width=1200):
    """
    Render a matplotlib figure as a fixed-width SVG.
    """
    buffer = io.BytesIO()

    fig.savefig(
        buffer,
        format="svg",
        bbox_inches="tight",
    )

    buffer.seek(0)

    svg_base64 = base64.b64encode(buffer.getvalue()).decode("utf-8")

    st.markdown(
        f"""
        <div style="width: {width}px; max-width: 100%; overflow-x: hidden;">
            <img
                src="data:image/svg+xml;base64,{svg_base64}"
                style="width: {width}px; max-width: 100%; height: auto;"
            />
        </div>
        """,
        unsafe_allow_html=True,
    )


@st.cache_data
def load_reference_data():
    """
    Load preprocessed HPA reference files.
    """
    if not CLUSTER_REFERENCE.exists():
        raise FileNotFoundError(
            f"Missing cluster reference file: {CLUSTER_REFERENCE}"
        )

    if not CELLTYPE_REFERENCE.exists():
        raise FileNotFoundError(
            f"Missing cell-type reference file: {CELLTYPE_REFERENCE}"
        )

    clusters = pd.read_parquet(CLUSTER_REFERENCE)
    hpa = pd.read_parquet(CELLTYPE_REFERENCE)

    return clusters, hpa


@st.cache_data
def load_available_tissues():
    """
    Load available tissues for display on the homepage and sidebar.
    """
    if not CLUSTER_REFERENCE.exists():
        raise FileNotFoundError(
            f"Missing cluster reference file: {CLUSTER_REFERENCE}"
        )

    clusters = pd.read_parquet(CLUSTER_REFERENCE)
    return list_tissues(clusters)


def read_gene_list(file_name, file_bytes):
    """
    Read uploaded candidate gene list.
    """
    buffer = io.BytesIO(file_bytes)
    file_name_lower = file_name.lower()

    if file_name_lower.endswith(".csv"):
        return pd.read_csv(buffer)

    return pd.read_csv(buffer, sep="\t")


def make_genes_passing_threshold_table(results):
    """
    Create a summary table listing each unique gene that passes the threshold.
    """
    expression_col = results.expression_col

    if results.filtered.empty:
        return pd.DataFrame(
            columns=[
                "Gene name",
                "Number of cell types",
                "Cell types passing threshold",
                f"Max {expression_col}",
                f"Mean {expression_col}",
            ]
        )

    genes_df = results.filtered.copy()

    summary_df = (
        genes_df.groupby("Gene name")
        .agg(
            **{
                "Number of cell types": ("Cell type", "nunique"),
                "Cell types passing threshold": (
                    "Cell type",
                    lambda cells: ", ".join(sorted(set(cells))),
                ),
                f"Max {expression_col}": (expression_col, "max"),
                f"Mean {expression_col}": (expression_col, "mean"),
            }
        )
        .reset_index()
        .sort_values(
            [f"Max {expression_col}", "Gene name"],
            ascending=[False, True],
        )
        .reset_index(drop=True)
    )

    summary_df[f"Max {expression_col}"] = summary_df[f"Max {expression_col}"].round(2)
    summary_df[f"Mean {expression_col}"] = summary_df[f"Mean {expression_col}"].round(2)

    return summary_df


def make_top_ncpm_plot(results, top_n=10, allowed_genes=None):
    """
    Plot the top cell-gene combinations by average nCPM.
    """
    if results.ncpm_df.empty:
        return None, pd.DataFrame()

    required_cols = ["Gene name", "Cell type", "nCPM"]

    missing_cols = [
        col for col in required_cols
        if col not in results.ncpm_df.columns
    ]

    if missing_cols:
        raise ValueError(f"nCPM table is missing columns: {missing_cols}")

    plot_df = results.ncpm_df[required_cols].copy()

    if allowed_genes is not None:
        plot_df = plot_df[plot_df["Gene name"].isin(allowed_genes)].copy()

    plot_df["nCPM"] = pd.to_numeric(plot_df["nCPM"], errors="coerce")
    plot_df = plot_df.dropna(subset=["Gene name", "Cell type", "nCPM"])

    if plot_df.empty:
        return None, pd.DataFrame()

    top_df = (
        plot_df.groupby(["Cell type", "Gene name"], as_index=False)["nCPM"]
        .mean()
        .rename(columns={"nCPM": "Average nCPM"})
        .sort_values("Average nCPM", ascending=False)
        .head(int(top_n))
        .reset_index(drop=True)
    )

    if top_df.empty:
        return None, top_df

    top_df["Cell type + gene"] = (
        top_df["Cell type"].astype(str)
        + " | "
        + top_df["Gene name"].astype(str)
    )

    fig, ax = plt.subplots(figsize=(18, 10))

    ax.bar(top_df["Cell type + gene"], top_df["Average nCPM"])

    ax.set_xlabel("Cell type and gene", fontsize=14)
    ax.set_ylabel("Average nCPM", fontsize=14)
    ax.set_title(
        f"Top {top_n} cell-gene combinations by average nCPM",
        fontsize=18,
        pad=18,
    )

    ax.tick_params(axis="x", rotation=65, labelsize=12)
    ax.tick_params(axis="y", labelsize=12)

    for label in ax.get_xticklabels():
        label.set_ha("right")

    ax.margins(x=0.005)

    fig.subplots_adjust(
        left=0.07,
        right=0.99,
        bottom=0.46,
        top=0.90,
    )

    return fig, top_df


# -----------------------------
# Sidebar
# -----------------------------

st.sidebar.header("1. Upload gene list")

use_example_gene_list = st.sidebar.checkbox(
    "Use example gene list",
    value=False,
)

gene_file = st.sidebar.file_uploader(
    "Upload candidate gene list TSV, TXT, or CSV",
    type=["tsv", "txt", "csv"],
    disabled=use_example_gene_list,
)

st.sidebar.header("2. Choose settings")

try:
    tissue_options = load_available_tissues()
except Exception:
    tissue_options = ["Immune cells"]

default_tissue_index = 0

if "Immune cells" in tissue_options:
    default_tissue_index = tissue_options.index("Immune cells")

selected_tissue = st.sidebar.selectbox(
    "Choose tissue",
    options=tissue_options,
    index=default_tissue_index,
)

threshold = st.sidebar.number_input(
    "Selected-cell expression threshold",
    min_value=0.0,
    value=1.0,
    step=0.5,
)

non_selected_threshold = st.sidebar.number_input(
    "Expression threshold for other cell types",
    min_value=0.0,
    value=float(threshold),
    step=0.5,
)

max_non_selected_cell_types = st.sidebar.number_input(
    "Maximum number of other cell types allowed above threshold",
    min_value=0,
    max_value=50,
    value=3,
    step=1,
)

use_selective_genes_for_plot = st.sidebar.checkbox(
    "Plot only selective genes",
    value=True,
)

top_n = st.sidebar.number_input(
    "Number of gene-cell type combinations to show",
    min_value=1,
    max_value=100,
    value=10,
    step=1,
)

run_button = st.sidebar.button("Run MendelCell analysis")


# -----------------------------
# Reference file status
# -----------------------------

with st.expander("Reference file status"):
    st.write("Expected reference files:")

    st.code(str(CLUSTER_REFERENCE))

    if CLUSTER_REFERENCE.exists():
        st.success(
            f"Cluster reference found "
            f"({CLUSTER_REFERENCE.stat().st_size / 1_000_000:.2f} MB)"
        )
    else:
        st.error("Cluster reference file is missing.")

    st.code(str(CELLTYPE_REFERENCE))

    if CELLTYPE_REFERENCE.exists():
        st.success(
            f"Cell-type reference found "
            f"({CELLTYPE_REFERENCE.stat().st_size / 1_000_000:.2f} MB)"
        )
    else:
        st.error("Cell-type reference file is missing.")

    st.code(str(EXAMPLE_GENE_LIST))

    if EXAMPLE_GENE_LIST.exists():
        st.success("Example gene list found.")
    else:
        st.warning("Example gene list is missing.")


# -----------------------------
# Stop if no gene file uploaded
# -----------------------------

if gene_file is None and not use_example_gene_list:
    st.info(
        "Upload a candidate gene list to begin, "
        "or select **Use example gene list** in the sidebar."
    )

    st.markdown(
        """
        Your gene list should contain a column named:

        ```text
        Gene Symbol
        ```

        Example:

        ```text
        Gene Symbol
        INS
        GCG
        PDX1
        CD3D
        PTPRC
        ```

        To test the app on Hugging Face without uploading your own file, select:

        ```text
        Use example gene list
        ```
        """
    )

    st.header("Available tissues")

    st.write(
        "Choose one of these tissue names in the sidebar. "
        "You can also choose **Immune cells** to analyze immune-related cell types."
    )

    try:
        available_tissues = load_available_tissues()

        tissue_df = pd.DataFrame(
            {
                "Tissue name": available_tissues
            }
        )

        with st.expander("Show available tissue names", expanded=True):
            show_dataframe_with_1_index(tissue_df, height=300, width=600)

    except Exception as e:
        st.warning(f"Could not load available tissue names: {e}")

    st.stop()


# -----------------------------
# Read gene list
# -----------------------------

try:
    if use_example_gene_list:
        if not EXAMPLE_GENE_LIST.exists():
            raise FileNotFoundError(
                f"Missing example gene list file: {EXAMPLE_GENE_LIST}"
            )

        gene_table = pd.read_csv(EXAMPLE_GENE_LIST, sep="\t")
        gene_file_name = EXAMPLE_GENE_LIST.name

    else:
        gene_table = read_gene_list(gene_file.name, gene_file.getvalue())
        gene_file_name = gene_file.name

except Exception as e:
    st.error(f"Could not read gene list file: {e}")
    st.stop()


# -----------------------------
# Validate gene list
# -----------------------------

if "Gene Symbol" not in gene_table.columns:
    st.error("Gene list must contain a column named 'Gene Symbol'.")
    st.write("Columns found:")
    st.write(list(gene_table.columns))
    st.stop()


st.success(
    f"Gene list loaded: {gene_file_name} "
    f"({gene_table.shape[0]:,} rows and {gene_table.shape[1]:,} columns)"
)

with st.expander("Preview gene list"):
    show_dataframe_with_1_index(gene_table.head(20), height=250, width=900)


# -----------------------------
# Wait for Run button
# -----------------------------

if not run_button:
    st.info("Choose a tissue and threshold, then click **Run MendelCell analysis**.")
    st.stop()


# -----------------------------
# Load reference files after Run
# -----------------------------

try:
    with st.spinner("Loading MendelCell reference files..."):
        clusters, hpa = load_reference_data()

    st.success(
        f"Reference loaded: {clusters.shape[0]:,} cluster rows and "
        f"{hpa.shape[0]:,} cell-type rows"
    )

except Exception as e:
    st.error("Could not load MendelCell reference files.")
    st.error(str(e))
    st.stop()


# -----------------------------
# Run MendelCell analysis
# -----------------------------

try:
    with st.spinner("Running MendelCell analysis..."):
        results = run_mendelcell(
            clusters=clusters,
            hpa=hpa,
            gene_table=gene_table,
            tissue=selected_tissue,
            threshold=threshold,
            non_selected_threshold=non_selected_threshold,
            max_non_selected_cell_types=max_non_selected_cell_types,
        )

except Exception as e:
    st.error(f"MendelCell analysis failed: {e}")
    st.stop()


# -----------------------------
# Display results
# -----------------------------

st.success("Analysis complete.")

selective_genes_df = results.selective_genes_df

col1, col2, col3, col4 = st.columns(4)

col1.metric("Input genes", len(results.gene_symbols))
col2.metric("Tissue-specific cell types", len(results.unique_cells))
col3.metric("Genes passing threshold", results.filtered["Gene name"].nunique())
col4.metric(
    "Selective genes",
    len(selective_genes_df["Gene name"].unique()) if not selective_genes_df.empty else 0,
)


st.header("Genes passing threshold")

genes_passing_threshold_df = make_genes_passing_threshold_table(results)

if genes_passing_threshold_df.empty:
    st.info("No genes passed the selected threshold.")
else:
    show_dataframe_with_1_index(
        genes_passing_threshold_df,
        height=350,
        width=1200,
    )


st.header("Selective genes")

st.write(
    "These genes pass the selected-cell expression threshold and are allowed "
    "to be above the other-cell threshold in only a limited number of other cell types."
)

if selective_genes_df.empty:
    st.info("No selective genes were found using the current thresholds.")
else:
    show_dataframe_with_1_index(
        selective_genes_df,
        height=350,
        width=1200,
    )


if use_selective_genes_for_plot:
    allowed_genes = set(selective_genes_df["Gene name"])
    plot_header = f"Top {top_n} selective cell-gene combinations by average nCPM"
else:
    allowed_genes = None
    plot_header = f"Top {top_n} cell-gene combinations by average nCPM"


st.header(plot_header)

try:
    top_ncpm_fig, top_ncpm_df = make_top_ncpm_plot(
        results,
        top_n=top_n,
        allowed_genes=allowed_genes,
    )

    if top_ncpm_fig is None or top_ncpm_df.empty:
        st.info(f"No nCPM values available for top-{top_n} plotting.")
    else:
        show_matplotlib_svg(top_ncpm_fig, width=1200)
        plt.close(top_ncpm_fig)

        show_dataframe_with_1_index(top_ncpm_df, height=400, width=1200)

except Exception as e:
    st.error(f"Could not create top-{top_n} nCPM plot.")
    st.stop()


# -----------------------------
# Download outputs
# -----------------------------

st.header("Download outputs")

unique_tsv = results.unique_to_tissue.to_csv(sep="\t", index=False)
filtered_tsv = results.filtered.to_csv(sep="\t", index=False)
ncpm_tsv = results.ncpm_df.to_csv(sep="\t", index=False)

if not genes_passing_threshold_df.empty:
    genes_passing_threshold_tsv = genes_passing_threshold_df.to_csv(
        sep="\t",
        index=False,
    )
else:
    genes_passing_threshold_tsv = ""

if not selective_genes_df.empty:
    selective_genes_tsv = selective_genes_df.to_csv(
        sep="\t",
        index=False,
    )
else:
    selective_genes_tsv = ""

if "top_ncpm_df" in locals() and not top_ncpm_df.empty:
    top_ncpm_tsv = top_ncpm_df.to_csv(sep="\t", index=False)
else:
    top_ncpm_tsv = ""

safe_tissue = safe_filename(results.selected_tissue)
pdf_filename = f"MendelCell_report_{safe_tissue}.pdf"


try:
    with tempfile.TemporaryDirectory() as tmpdir:
        pdf_path = Path(tmpdir) / pdf_filename
        create_pdf_report(results, pdf_path)

        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()

    st.download_button(
        label="Download PDF report",
        data=pdf_bytes,
        file_name=pdf_filename,
        mime="application/pdf",
    )

except Exception as e:
    st.error("Could not create PDF report.")
    st.error(str(e))


if genes_passing_threshold_tsv:
    st.download_button(
        label="Download genes passing threshold TSV",
        data=genes_passing_threshold_tsv,
        file_name="genes_passing_threshold.tsv",
        mime="text/tab-separated-values",
    )

if selective_genes_tsv:
    st.download_button(
        label="Download selective genes TSV",
        data=selective_genes_tsv,
        file_name="selective_genes.tsv",
        mime="text/tab-separated-values",
    )

st.download_button(
    label="Download unique cell types TSV",
    data=unique_tsv,
    file_name="unique_cell_types.tsv",
    mime="text/tab-separated-values",
)

st.download_button(
    label="Download filtered candidate genes TSV",
    data=filtered_tsv,
    file_name="filtered_candidate_genes.tsv",
    mime="text/tab-separated-values",
)

st.download_button(
    label="Download nCPM table TSV",
    data=ncpm_tsv,
    file_name="candidate_gene_ncpm_by_cell_type.tsv",
    mime="text/tab-separated-values",
)

if top_ncpm_tsv:
    st.download_button(
        label=f"Download top {top_n} average nCPM TSV",
        data=top_ncpm_tsv,
        file_name=f"top_{top_n}_cell_gene_average_ncpm.tsv",
        mime="text/tab-separated-values",
    )