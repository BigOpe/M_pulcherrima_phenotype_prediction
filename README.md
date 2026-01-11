Structured Phenotypic Architecture and Functional Genomic Signatures in Metschnikowia pulcherrima
This repository contains the complete analysis pipeline used in the study:

Integrative Analysis of Phenotypic Traits and Genomic Features in Metschnikowia pulcherrima Using Machine Learning

The project integrates:

Phenotypic trait measurements across 15 M. pulcherrima strains

Genome annotation and functional categorization

Trait–trait predictive modelling (PLS regression, Random Forest)

Correlation networks and PCA

Exploratory dna2vec sequence embeddings

A past–present–future phenotypic projection model (2015–2035)

The overarching goal is to characterize structured phenotypic relationships in M. pulcherrima and provide a reproducible framework for exploratory genotype–phenotype interpretation.

📁 Repository Structure
Code
analysis/
    Core Python scripts for phenotypic modelling, PCA, projections, and figure generation.

dna2vec-legacy/
    Pretrained dna2vec embeddings used for exploratory sequence-level PCA.

data/  (if added)
    Phenotypic tables, contig lists, or annotation summaries.

environment.yml
    Conda environment file for full reproducibility.

.gitignore
    Ensures virtual-environment files and binaries are excluded.
⚙️ Environment Setup
Create the software environment used in this project:

bash
conda env create -f environment.yml
conda activate mpulch_ml
This installs all required packages (NumPy, pandas, scikit‑learn, seaborn, matplotlib, networkx, etc.).

▶️ How to Reproduce the Analysis
1. Create the Conda environment
bash
conda env create -f environment.yml
conda activate mpulcherrima_env
2. Run the Analysis Pipeline
Each script in analysis/ corresponds to a specific component of the study:

start.py — main pipeline entry point

animated_pca_.py — PCA projection and animation

gene.py / prok.py — genome annotation parsing

rel.py — regression modelling (PLS, Random Forest)

contig.py — assembly parsing

checking.py — QC and preprocessing

che.py, chh.py, etc. — supporting analysis modules

You may run scripts individually, or execute the full workflow using:

bash
bash run_all.sh
This script performs:

Genome annotation and functional categorization

Phenotypic data standardization

Exploratory gene–trait mapping

Correlation analysis and network construction

Principal Component Analysis (PCA)

Machine‑learning regression (PLS, Random Forest)

Phenotypic divergence projection (2015–2035)

📊 Data Availability
All phenotypic data, annotation summaries, and analysis scripts required to reproduce the figures and results are included in this repository.
The full computational workflow is reproducible using the provided environment.yml.

🧪 Reproducibility Notes
All analyses were performed in Python using a fully specified Conda environment.

No virtual‑environment files or binaries are included in the repository.

All scripts are self‑contained and can be run independently.

dna2vec embeddings are used only for exploratory PCA, not as predictors in ML models.

All random processes use fixed random seeds.

Preprocessing steps are confined within cross‑validation loops to prevent data leakage.

Genomic analyses are exploratory and not used as predictors in machine‑learning models.

📬 Contact
For questions or collaboration inquiries:

Opeyemi Adesuyi  
Bioinformatics and Molecular Biology Researcher
Email: ofadesuyi@unite.it
