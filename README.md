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

MendelCell

MendelCell is an app that helps prioritize candidate genes by looking at their expression across tissue-specific and cell-type-specific single-cell reference data.

The app was built to support exploratory analysis of candidate genes from rare disease, monogenic disease, and immune-related genetics projects.

What MendelCell Does

MendelCell allows users to upload a candidate gene list and identify which genes are expressed in relevant cell types.

The app can:

Upload a candidate gene list as a TSV, TXT, or CSV file
Analyze a selected tissue, such as Pancreas
Analyze the special pseudo-tissue Immune cells
Apply an expression threshold
Show candidate genes found in each cell type
Generate summary tables
Generate nCPM expression plots
Show the top 10 gene-cell type combinations by average nCPM
Download TSV results and a PDF report
Input File Format

The uploaded file must contain a column named:

Gene Symbol

Example:

Gene Symbol
INS
GCG
PDX1
CD3D
PTPRC
IL2RA
CTLA4
Tissue Options

You can enter tissue names such as:

Pancreas
Liver
Kidney
Lung

You can also enter:

Immune cells

or:

immune

This runs the analysis on immune-related cell types across the reference data.

Output

MendelCell generates:

Cell types associated with the selected tissue
Candidate gene counts per cell type
Candidate genes detected in each cell type
Filtered gene-cell type expression results
Mean nCPM expression values
A top 10 average nCPM plot
Downloadable TSV files
A downloadable PDF report