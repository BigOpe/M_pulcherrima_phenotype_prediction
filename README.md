# Metschnikowia pulcherrima phenotype prediction

This repository contains code associated with a machine-learning analysis
of phenotypic traits and genomic features in *Metschnikowia pulcherrima* strains.

# Structured Phenotypic Architecture and Functional Genomic Signatures in *Metschnikowia pulcherrima*

This repository contains the full analysis pipeline used in the study:

The project integrates:
- Phenotypic trait measurements across 15 strains
- Genome annotation and functional categorization
- Trait–trait predictive modelling (PLS, Random Forest)
- Correlation networks and PCA
- Exploratory dna2vec sequence embeddings
- A past–present–future projection model (2015–2035)

The goal is to characterize structured phenotypic relationships in *M. pulcherrima* and provide a reproducible framework for exploratory genotype–phenotype interpretation.

---

## 📁 Repository Structure

analysis/
Core Python scripts for phenotypic modelling, PCA, projections, and figure generation.

dna2vec-legacy/
Pretrained dna2vec embeddings used for exploratory sequence-level PCA.

data/ (if added)
Phenotypic tables, contig lists, or annotation summaries.

environment.yml
Conda environment file for full reproducibility.

## Environment setup

To create the software environment used in this project:

```bash
conda env create -f environment.yml
conda activate mpulch_ml
.gitignore
Ensures virtual-environment files and binaries are excluded.


---

## ▶️ How to Reproduce the Analysis

### **1. Create the Conda environment**

conda env create -f environment.yml
conda activate mpulcherrima_env


This installs all required packages (NumPy, pandas, scikit‑learn, seaborn, matplotlib, networkx, etc.).

---

### **2. Run the analysis scripts**

Each script in `analysis/` corresponds to a specific part of the study:

- `start.py` — main pipeline entry point  
- `animated_pca_.py` — PCA projection and animation  
- `gene.py` / `prok.py` — genome annotation parsing  
- `rel.py` — regression modelling (PLS, Random Forest)  
- `contig.py` — assembly parsing  
- `checking.py` — QC and preprocessing  
- `che.py`, `chh.py`, etc. — supporting analysis modules  

Run any script individually:

python analysis/start.py



or execute scripts step‑by‑step depending on the analysis you want to reproduce.

---

## 📊 Data Availability

All phenotypic data, annotation summaries, and analysis scripts required to reproduce the figures and results are included in this repository.

The full computational workflow is reproducible using the provided `environment.yml`.

---


## 🧪 Reproducibility Notes

- All analyses were performed in Python using a fully specified Conda environment.  
- No virtual‑environment files or binaries are included in the repository.  
- All scripts are self‑contained and can be run independently.  
- dna2vec embeddings are included only for exploratory PCA and are not used as predictors in ML models.

---

## 📬 Contact

For questions or collaboration inquiries, please contact:

**Opeyemi Adesuyi**  
Bioinformatics and Molecular Biology Researcher  
Email: [ofadesuyi@unite.it]



