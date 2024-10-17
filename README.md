# AsaruSim paper scripts
The current repository contains all the scripts needed to reproduce the results published in the paper:
> Ali Hamraoui, Laurent Jourdren and Morgane Thomas-Chollier. AsaruSim: a single-cell and spatial RNA-Seq Nanopore long-reads simulation workflow.
> bioRxiv 2024.09.20.613625; doi: [https://doi.org/10.1101/2024.09.20.613625](https://doi.org/10.1101/2024.09.20.613625)

## Structure

The repository is organized into three main directories, each pertaining to a different aspect of RNA-Seq simulation and analysis:

### 1. Sequence Identity Estimation

This directory contains scripts and notebooks for estimating sequence identity between simulated and real reads.

- `identity.html` - HTML output of the Jupyter Notebook analysis.
- `identity.ipynb` - Jupyter Notebook containing the analysis workflows.
- `optimization.py` - Python script for optimizing simulation parameters.
- `seqAligner.py` - Python script for sequence alignment processes.

### 2. Read Length Distribution

Tools and visualizations for estimating and comparing the read length distributions of simulated and real RNA-Seq reads.

- `lognorm.pdf` - PDF file showing the log-normal distribution fit of read lengths.
- `read_length_estimation.html` - HTML output of the read length estimation Jupyter Notebook.
- `read_length_estimation.ipynb` - Jupyter Notebook detailing the estimation of read lengths.
- `real_vs_annotation-4.pdf` - PDF comparison between real reads and annotations.
- `real_vs_simulation.pdf` - PDF showing comparison of read lengths between real and simulated data.

### 3. Real vs. Simulated Expression

Analysis of gene expression differences between real and simulated datasets.

- `Real_PBMC.pdf` - Visualization of real PBMC dataset expression levels.
- `Simulated_PBMC.pdf` - Visualization of simulated PBMC dataset expression levels.
- `TSNE_milisi.pdf` - t-SNE plot comparing real and simulated datasets.
- `UMAP_comparaison.Rmd` - R Markdown document for UMAP analysis comparing real and simulated expression.
- `UMAP_comparaison.html` - HTML output of the UMAP comparison.
- `imports.R` - R script containing libraries necessary for the analysis.

## Getting Started

To get started with AsaruSim, clone this repository and explore the notebooks and scripts provided. Each directory includes detailed notebooks and scripts necessary for performing specific types of analysis. Ensure you have the required dependencies installed, which can be identified in the respective scripts and notebooks.

## Requirements

Make sure to install the necessary Python and R libraries before running the scripts. Dependencies include but are not limited to libraries like `numpy`, `matplotlib`, `seaborn` for Python, and `ggplot2`, `dplyr`, `cowplot` for R. You can install Python libraries via `pip` and R packages via `install.packages()`.

## Contributing

Contributions to AsaruSim are welcome! Please feel free to fork the repository, make your changes, and submit a pull request. If you have any questions or need further information, do not hesitate to raise an issue or contact the repository maintainers.

Thank you for visiting the AsaruSim repository, and we look forward to your contributions and feedback to enhance this simulation workflow!
