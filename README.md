# Glioblastoma ERβ Isoform Explorer

An interactive R/Shiny dashboard designed to visualize transcriptional changes in U87 glioblastoma cells. This tool compares gene expression across Control, ERβ-Knockout, and the reintroduction of ERβ1 and ERβ5 isoforms.

## Biological Context
This app utilizes RNA-seq data (TPM) from study **GSE104296**. It explores how specific Estrogen Receptor beta (ERβ) isoforms modulate unique pathways such as NF-κB, Jak/STAT, and mTOR signaling.

##  Features
- **Dynamic Search:** Explore 39,000+ genes using NCBI IDs or Gene Symbols.
- **Statistical Inference:** Automated p-value calculation comparing experimental groups to Control.
- **Data Integration:** Integrated raw sequencing counts with Bioconductor's `org.Hs.eg.db`.
- **Interactive UI:** Built with `shiny`, `ggplot2`, and `DT`.

## How to Run Locally
1. Clone this repository.
2. Ensure you have R/RStudio installed.
3. Install dependencies:
   `install.packages(c("shiny", "ggplot2", "dplyr", "DT", "ggpubr", "readr"))`
4. Run `shiny::runApp()`
