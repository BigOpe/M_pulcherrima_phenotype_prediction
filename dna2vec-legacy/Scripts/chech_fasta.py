import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_predict, KFold, StratifiedKFold
from sklearn.metrics import r2_score, accuracy_score, ConfusionMatrixDisplay
from adjustText import adjust_text
from Bio import SeqIO
import itertools
from collections import Counter
import umap

# ------------------------------
# File paths
# ------------------------------
phen_file   = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
fasta_file  = r"C:\Users\Giorgia\dna\dna2vec-legacy\gen.fasta"
map_file    = r"C:\Users\Giorgia\dna\predicted_accessions.csv"

# ------------------------------
# Utility Functions
# ------------------------------

def clean_columns(df):
    df.columns = (
        df.columns.str.replace(r"\(.*?\)", "", regex=True)
                  .str.replace("% v/v", "", regex=False)
                  .str.replace(" ", "_")
                  .str.lower()
                  .str.strip("_")
    )
    return df

def map_categorical(df, col, mapping):
    return df[col].map(mapping)

def get_kmer_freqs(seq, k=4):
    counts = Counter(seq[i:i+k] for i in range(len(seq)-k+1) if "N" not in seq[i:i+k])
    total = sum(counts.values())
    return [counts.get(kmer, 0)/total if total>0 else 0 for kmer in kmers]

def plot_2D(X_2D, labels=None, title="", palette=None, size=30):
    plt.figure(figsize=(8,6))
    if labels is not None and palette is not None:
        unique_labels = sorted(set(labels))
        for i, lbl in enumerate(unique_labels):
            idx = [j for j, l in enumerate(labels) if l == lbl]
            plt.scatter(X_2D[idx,0], X_2D[idx,1], s=size, color=palette[i], label=lbl, alpha=0.8)
        plt.legend(bbox_to_anchor=(1.05,1), loc="upper left", fontsize=8)
    else:
        plt.scatter(X_2D[:,0], X_2D[:,1], s=size, color='steelblue', alpha=0.7)
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")
    plt.title(title)
    plt.tight_layout()
    plt.show()

# ------------------------------
# Load and clean phenotype data
# ------------------------------
df_pheno = pd.read_csv(phen_file)
df_pheno = clean_columns(df_pheno)
df_pheno = df_pheno.dropna(subset=["strain"])  # Ensure strain exists

# ------------------------------
# Map categorical traits to numeric
# ------------------------------
categorical_mappings = {
    "killer_activity": {"Weak":0.2,"Moderate":0.6,"Strong":1.0},
    "pigment_intensity": {"Moderate":0.6,"High":0.9},
    "ethanol_tolerance": lambda x: float(x)/10 if pd.notnull(x) else np.nan,
    "copper_resistance": lambda x: float(x)/4 if pd.notnull(x) else np.nan
}

for col, mapping in categorical_mappings.items():
    if callable(mapping):
        df_pheno[f"{col}_num"] = df_pheno[col].apply(mapping)
    else:
        df_pheno[f"{col}_num"] = df_pheno[col].map(mapping)

numeric_cols = [c for c in df_pheno.columns if c.endswith("_num")]

# ------------------------------
# PCA on numeric traits
# ------------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df_pheno[numeric_cols])
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X_scaled)
df_pheno[["PC1","PC2"]] = X_pca
print("Explained variance ratio (PCA):", pca.explained_variance_ratio_)

plot_2D(X_pca, labels=df_pheno["strain"], title="PCA of Numeric Traits",
        palette=sns.color_palette("tab20", len(df_pheno["strain"].unique())))

# ------------------------------
# Regression/Classification per trait
# ------------------------------
for target in numeric_cols:
    X = df_pheno[[c for c in numeric_cols if c != target]].values
    y = df_pheno[target].values

    model = Ridge(alpha=1.0)
    cv = KFold(n_splits=min(5, len(df_pheno)))
    y_pred = cross_val_predict(model, X, y.astype(float), cv=cv)
    r2 = r2_score(y, y_pred)

    plt.figure(figsize=(8,6))
    plt.scatter(y, y_pred, s=40, color='steelblue', alpha=0.7)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], "r--")
    plt.xlabel("Observed")
    plt.ylabel("Predicted")
    plt.title(f"{target} regression\nR²={r2:.2f}")
    plt.tight_layout()
    plt.show()

# ------------------------------
# Load FASTA & build k-mer matrix
# ------------------------------
kmers = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
sequences, accessions = [], []

for rec in SeqIO.parse(fasta_file, "fasta"):
    accessions.append(rec.id.split()[0].split(".")[0].upper())
    sequences.append(str(rec.seq).upper())

X_kmer = np.array([get_kmer_freqs(seq) for seq in sequences])
X_kmer_scaled = StandardScaler().fit_transform(X_kmer)

# PCA & UMAP
X_pca_seq = PCA(n_components=2).fit_transform(X_kmer_scaled)
X_umap_seq = umap.UMAP(random_state=42).fit_transform(X_kmer_scaled)

plot_2D(X_pca_seq, title="PCA of Accessions (k-mer features)")
plot_2D(X_umap_seq, title="UMAP of Accessions (k-mer features)", color='forestgreen')

# ------------------------------
# Map accessions to strains
# ------------------------------
df_map = pd.read_csv(map_file)
df_map = clean_columns(df_map)
df_map["accession"] = df_map["accession"].astype(str).str.upper().str.split(".").str[0]
map_dict = dict(zip(df_map["accession"], df_map["strain"]))

mapped_strains = [map_dict.get(acc, "UNMAPPED") for acc in accessions]
unmapped = mapped_strains.count("UNMAPPED")
print(f"Mapped: {len(accessions)-unmapped}, Unmapped: {unmapped}")
