# BNMF Genomic Clustering Pipeline

A professional Python/R hybrid pipeline designed for processing GWAS summary statistics prodused by postgwas harmonisation and performing Bayesian Non-negative Matrix Factorization (bNMF) clustering. This tool implements the partitioning methodology described in the [gwas-partitioning/bnmf-clustering](https://github.com/gwas-partitioning/bnmf-clustering) repository.

## 🏗 Project Structure & Workflow

The pipeline automates the transition from raw GWAS summary statistics to biological clusters via three integrated steps:

* **Step 1: Main GWAS Processing (Python/Polars)** Extracts independent lead variants from the primary GWAS summary statistics. It utilizes high-performance LD clumping to identify the most significant, non-redundant genetic signals within defined genomic windows.
* **Step 2: Trait Matrix Construction (Python/VCF)** Generates synchronized Z-score and effective sample size () matrices. This step queries secondary trait GWAS files (VCF/CSV) at the specific lead SNP positions identified in Step 1 to build the input matrix for factorization.
* **Step 3: bNMF Analysis (R/rpy2)** Executes the core clustering algorithm. It performs stability analysis across multiple repetitions (default: 100 reps) to determine the optimal number of clusters () and generates automated HTML reports summarizing the biological traits associated with each cluster.

## 🚀 Installation

### 1. Environment Setup

It is highly recommended to use the provided `conda_env.yaml` to ensure all C-libraries (like `bcftools`) and R dependencies are correctly linked to your Python environment.

```bash
# Create the environment from the project root
conda env create -f bnmf_pkg/conda_env.yaml

# Activate the environment
conda activate bnmf_pkg_env

```

### 2. Package Installation

Install the project as an editable Python package to register the `bnmf-pipeline` command globally.

```bash
# Navigate to the root directory (containing pyproject.toml)
cd /Users/JJOHN41/Documents/developing_software/bnmf_pkg

# Install in editable mode
pip install -e .

```

## 💻 Usage

The pipeline is controlled via a single command-line interface with modern, color-coded formatting.

### Commands

* **Full Pipeline:** `bnmf-pipeline --run_traits --run_main --run_bnmf`
* **Custom Repetitions:** `bnmf-pipeline --run_bnmf --n_reps 250 --main_gwas_id "SCZ_2026"`
* **View Documentation:** `bnmf-pipeline --help`

## 📋 Requirements

* **Python:** 3.12+ (`polars`, `rpy2`, `rich`, `rich-argparse`)
* **R:** 4.4+ (`tidyverse`, `data.table`, `GenomicRanges`, `rmarkdown`)
* **System Tools:** `bcftools` must be accessible in your system `$PATH`.

## 📄 References

* **Methodology:** [gwas-partitioning/bnmf-clustering](https://github.com/gwas-partitioning/bnmf-clustering)

---

Would you like me to help you verify that the `pyproject.toml` is correctly configured to support the `bnmf-pipeline` command once the package is installed?