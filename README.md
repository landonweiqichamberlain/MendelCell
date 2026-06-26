---
title: MendelCell
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 8501
pinned: false
---
---

title: MendelCell
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 8501
pinned: false
-------------

# MendelCell

**MendelCell** is a Streamlit-based bioinformatics tool for prioritizing candidate genes by tissue-specific and cell-type-specific expression patterns.

The app uses a preprocessed Human Protein Atlas single-cell reference to help identify which candidate genes are expressed in relevant cell types. MendelCell was designed for exploratory analysis of rare disease and monogenic disease candidate genes, with an emphasis on immune and pancreatic cell-type biology.

## Live App

The app can be deployed using Hugging Face Spaces or run locally with Streamlit.

Example Hugging Face Space URL:

```text
https://landonchamberlain-mendelcell.hf.space
```

## Features

MendelCell allows users to:

* Upload a candidate gene list as a TSV, TXT, or CSV file
* Select a tissue or pseudo-tissue, such as `Immune cells`
* Apply an expression threshold
* Identify tissue-specific or immune-relevant cell types
* View candidate genes expressed in each cell type
* Generate summary tables
* Generate plots of candidate gene expression
* Plot the top 10 gene-cell type combinations by average nCPM
* Download TSV outputs
* Download a PDF report

## Input Format

The uploaded gene list must contain a column named:

```text
Gene Symbol
```

Example TSV file:

```text
Gene Symbol
INS
GCG
PDX1
CD3D
PTPRC
IL2RA
CTLA4
```

CSV files are also supported, as long as they contain a `Gene Symbol` column.

## Tissue Selection

The app can analyze real tissues from the preprocessed reference, such as:

```text
Pancreas
Lung
Liver
Kidney
```

MendelCell also includes a special pseudo-tissue:

```text
Immune cells
```

Typing `Immune cells` or `immune` will analyze immune-related cell types across the reference data.

## Output Tables

MendelCell generates several outputs:

### Cell types unique to selected tissue

Lists cell types associated with the selected tissue or pseudo-tissue.

### Candidate gene count per cell type

Shows how many uploaded candidate genes are expressed in each cell type.

### Candidate genes found in each cell type

Lists candidate genes detected above the selected expression threshold.

### Filtered candidate genes

Shows gene-cell type pairs that pass the expression threshold.

### Mean nCPM values

Shows nCPM expression values for candidate genes in relevant cell types.

### Top 10 average nCPM plot

Ranks the top 10 gene-cell type combinations by average nCPM value, with the highest values shown first.

## Local Installation

Clone the repository:

```bash
git clone https://github.com/LandonChamberlain/MendelCell.git
cd MendelCell
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e .
```

Run the app:

```bash
python -m streamlit run mendelcell/app.py
```

## Requirements

The app requires:

```text
streamlit
pandas
matplotlib
pyarrow
```

The `requirements.txt` file should include:

```text
-e .
streamlit
pandas
matplotlib
pyarrow
```

## Project Structure

```text
MendelCell/
├── .github/
│   └── workflows/
│       └── sync-to-huggingface.yml
├── .streamlit/
│   └── config.toml
├── data/
│   ├── mendelcell_clusters_reference.parquet
│   └── mendelcell_celltype_reference.parquet
├── mendelcell/
│   ├── __init__.py
│   ├── app.py
│   ├── analysis.py
│   ├── io.py
│   └── report.py
├── scripts/
│   └── build_reference.py
├── Dockerfile
├── README.md
├── requirements.txt
├── pyproject.toml
└── .gitignore
```

## Reference Data

MendelCell uses preprocessed Human Protein Atlas single-cell expression reference files:

```text
data/mendelcell_clusters_reference.parquet
data/mendelcell_celltype_reference.parquet
```

The raw Human Protein Atlas TSV files are not included in the repository because they are large. Instead, MendelCell uses smaller preprocessed Parquet reference files for faster loading and deployment.

Raw files that should not be committed:

```text
rna_single_cell_cluster.tsv
rna_single_cell_type.tsv
rna_single_cell_cluster.tsv.zip
rna_single_cell_type.tsv.zip
```

## Hugging Face Deployment

This repository can be deployed as a Hugging Face Docker Space.

The `Dockerfile` runs:

```bash
python -m streamlit run mendelcell/app.py
```

The top of this README contains the Hugging Face Space metadata:

```yaml
---
title: MendelCell
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
app_port: 8501
pinned: false
---
```

## GitHub to Hugging Face Sync

This repository can be synced automatically to Hugging Face Spaces using GitHub Actions.

The workflow file is located at:

```text
.github/workflows/sync-to-huggingface.yml
```

After each push to the `main` branch, GitHub Actions can sync the latest code to the Hugging Face Space.

## Files Not Tracked by Git

The `.gitignore` file excludes local environments, Python cache files, generated metadata, large raw data files, and generated reports.

Examples:

```text
.venv/
__pycache__/
*.pyc
*.egg-info/
mendelcell.egg-info/
.DS_Store
*.pdf
mendelcell_output/
rna_single_cell_cluster.tsv
rna_single_cell_type.tsv
```

## Important Notes

MendelCell is intended for research and exploratory analysis only.

It is not a diagnostic tool and should not be used to make clinical decisions without appropriate validation and expert review.

## License

MIT License


## Citation

If you use MendelCell in research, please cite the GitHub repository and any relevant Human Protein Atlas data sources used to build the reference files.

## Author

Landon Chamberlain

Dublin High School
UCSF Diabetes Center research volunteer
