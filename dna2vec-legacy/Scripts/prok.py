# =========================================================
# Prokka Gene–Trait Functional Enrichment Visualization
# =========================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact

# -----------------------------
# 1. INPUT FILE
# -----------------------------
# This file should already exist from your previous steps
# Columns REQUIRED:
# locus_tag, product, strain,
# Ethanol_Tolerance, Killer_Activity, Pigment_Intensity,
# Copper_Resistance, Sugar_Utilization, H2S, Protease,
# Lipase, Growth_Rate, Glucosidase, SO2_Resistance,
# Esterase, Gluconic_Acid_Utilization

INPUT_FILE = r"C:\Users\Giorgia\dna\all_strains_gene_trait_matrix.csv"

df = pd.read_csv(INPUT_FILE)

# -----------------------------
# 2. DEFINE FUNCTIONAL GROUPS
# -----------------------------
# Map Prokka product descriptions → biological functions

FUNCTION_MAP = {
    "dehydrogenase": "Redox / Energy metabolism",
    "oxidoreductase": "Redox / Energy metabolism",
    "reductase": "Redox / Energy metabolism",
    "alcohol": "Ethanol metabolism",
    "transport": "Transporters",
    "permease": "Transporters",
    "ABC": "Transporters",
    "cell wall": "Cell wall organization",
    "mannoprotein": "Cell wall organization",
    "glucan": "Cell wall organization",
    "stress": "Stress response",
    "heat shock": "Stress response",
    "Hsp": "Stress response",
    "lipase": "Lipid metabolism",
    "esterase": "Lipid metabolism",
    "protease": "Protein degradation",
    "ubiquitin": "Protein turnover",
    "DNA": "Indeterminate / housekeeping",
    "RNA": "Indeterminate / housekeeping"
}

def assign_function(product):
    product = str(product).lower()
    for key, value in FUNCTION_MAP.items():
        if key.lower() in product:
            return value
    return "Hypothetical / unknown"

df["functional_category"] = df["product"].apply(assign_function)

# -----------------------------
# 3. SELECT TRAIT TO ANALYZE
# -----------------------------
TRAIT = "Ethanol_Tolerance"   # change freely

trait_genes = df[df[TRAIT] == 1]
background_genes = df[df[TRAIT] == 0]

# -----------------------------
# 4. ENRICHMENT ANALYSIS
# -----------------------------
results = []

for func in df["functional_category"].unique():

    a = len(trait_genes[trait_genes["functional_category"] == func])
    b = len(trait_genes) - a
    c = len(background_genes[background_genes["functional_category"] == func])
    d = len(background_genes) - c

    if a + c == 0:
        continue

    _, p = fisher_exact([[a, b], [c, d]], alternative="greater")

    results.append({
        "functional_category": func,
        "trait_genes": a,
        "background_genes": c,
        "p_value": p
    })

enrichment_df = pd.DataFrame(results)
enrichment_df["-log10(FDR)"] = -np.log10(enrichment_df["p_value"] + 1e-10)

# Keep meaningful categories only
enrichment_df = enrichment_df[enrichment_df["trait_genes"] > 3]
enrichment_df = enrichment_df.sort_values("trait_genes", ascending=True)

# -----------------------------
# 5. PANEL A — BAR PLOT
# -----------------------------
plt.figure(figsize=(7, 4))
plt.barh(
    enrichment_df["functional_category"],
    enrichment_df["trait_genes"],
    color="black"
)

plt.xlabel("Number of trait-associated genes")
plt.title(f"Top functional categories associated with {TRAIT}")
plt.tight_layout()
plt.savefig("Figure_Panel_A_Functional_Importance.png", dpi=300)
plt.show()

# -----------------------------
# 6. PANEL B — GO-STYLE DOT PLOT
# -----------------------------
plt.figure(figsize=(7, 6))

sns.scatterplot(
    data=enrichment_df,
    x="trait_genes",
    y="functional_category",
    size="trait_genes",
    hue="-log10(FDR)",
    palette="Reds",
    sizes=(40, 350),
    edgecolor="black"
)

plt.xlabel("Number of genes")
plt.ylabel("Biological process")
plt.title(f"Functional enrichment of {TRAIT}-associated genes")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig("Figure_Panel_B_Functional_Enrichment.png", dpi=300)
plt.show()

# -----------------------------
# 7. SAVE TABLE FOR PAPER
# -----------------------------
enrichment_df.to_csv(
    f"{TRAIT}_functional_enrichment_table.csv",
    index=False
)

print("✔ Analysis complete")
print("✔ Figures saved")
print("✔ Enrichment table exported")

import pandas as pd

# Load combined gene–trait / Prokka annotation table
df = pd.read_csv(r"C:\Users\Giorgia\dna\all_strains_gene_trait_matrix.csv")


# 2. Prefix → Strain mapping
# -----------------------------
prefix_map = {
    "VFXK": "FL01",
    "CP0344": "APC 1.2",
    "JACBPP": "KIOM G15050",
    "ANFW": "277",
    "JAQFWG": "DBT012",
    "JAJMHS": "NRRL Y-48711",
    "MTJM": "AP47",
    "JAKTYM": "NRRL Y-48710",
    "JAJMCQ": "CBS 15341",
    "JAJMIJ": "NRRL Y-7111",
    "JAJMIQ": "CBS 10809T",
    "QBLL": "UCD127",
    "JAKTYT": "NRRL Y-48712",
    "NISE": "BATH1",
    "JAKTUL": "NRRL Y-27328"
}

# Extract prefix from contig IDs (letters before digits)
df["prefix"] = df["strain"].str.extract(r"^([A-Z]+)")
df["strain_name"] = df["prefix"].map(prefix_map)

# -----------------------------
# 3. Keyword → Functional category mapping
# -----------------------------
category_map = {
    "dehydrogenase": "Redox / Energy metabolism",
    "oxidoreductase": "Redox / Energy metabolism",
    "reductase": "Redox / Energy metabolism",
    "alcohol": "Ethanol metabolism",
    "transport": "Transporters",
    "permease": "Transporters",
    "ABC": "Transporters",
    "cell wall": "Cell wall organization",
    "mannoprotein": "Cell wall organization",
    "glucan": "Cell wall organization",
    "stress": "Stress response",
    "heat shock": "Stress response",
    "Hsp": "Stress response",
    "lipase": "Lipid metabolism",
    "esterase": "Lipid metabolism",
    "protease": "Protein degradation",
    "ubiquitin": "Protein turnover",
    "DNA": "Indeterminate / housekeeping",
    "RNA": "Indeterminate / housekeeping"
}

def classify(product):
    product = str(product).lower()
    for keyword, category in category_map.items():
        if keyword.lower() in product:
            return category
    return "Other"

df["functional_category"] = df["product"].apply(classify)

# -----------------------------
# 4. Build summary table
# -----------------------------
# Number of genes per category
gene_counts = df.groupby("functional_category")["locus_tag"].count().reset_index()
gene_counts.columns = ["Functional category", "Number of genes"]

# Number of strains per category
strain_counts = df.groupby("functional_category")["strain_name"].nunique().reset_index()
strain_counts.columns = ["Functional category", "Number of strains"]

# Merge
summary = pd.merge(gene_counts, strain_counts, on="Functional category")

# -----------------------------
# 5. Save and display
# -----------------------------
summary.to_csv("functional_category_summary.csv", index=False)
print(summary)

import pandas as pd

# -----------------------------
# 1. Load annotation CSV
# -----------------------------
# Replace with the path to your file
df = pd.read_csv(r"C:\Users\Giorgia\dna\all_strains_gene_trait_matrix.csv")

# -----------------------------
# 2. Prefix → Strain mapping
# -----------------------------
prefix_map = {
    "VFXK": "FL01",
    "CP0344": "APC 1.2",
    "JACBPP": "KIOM G15050",
    "ANFW": "277",
    "JAQFWG": "DBT012",
    "JAJMHS": "NRRL Y-48711",
    "MTJM": "AP47",
    "JAKTYM": "NRRL Y-48710",
    "JAJMCQ": "CBS 15341",
    "JAJMIJ": "NRRL Y-7111",
    "JAJMIQ": "CBS 10809T",
    "QBLL": "UCD127",
    "JAKTYT": "NRRL Y-48712",
    "NISE": "BATH1",
    "JAKTUL": "NRRL Y-27328"
}

# Extract prefix from contig IDs (letters before digits)
df["prefix"] = df["strain"].str.extract(r"^([A-Z]+)")
df["strain_name"] = df["prefix"].map(prefix_map)

# -----------------------------
# 3. Keyword → Functional category mapping
# -----------------------------
category_map = {
    "dehydrogenase": "Redox / Energy metabolism",
    "oxidoreductase": "Redox / Energy metabolism",
    "reductase": "Redox / Energy metabolism",
    "alcohol": "Ethanol metabolism",
    "transport": "Transporters",
    "permease": "Transporters",
    "ABC": "Transporters",
    "cell wall": "Cell wall organization",
    "mannoprotein": "Cell wall organization",
    "glucan": "Cell wall organization",
    "stress": "Stress response",
    "heat shock": "Stress response",
    "Hsp": "Stress response",
    "lipase": "Lipid metabolism",
    "esterase": "Lipid metabolism",
    "protease": "Protein degradation",
    "ubiquitin": "Protein turnover",
    "DNA": "Indeterminate / housekeeping",
    "RNA": "Indeterminate / housekeeping"
}

def classify(product):
    product = str(product).lower()
    for keyword, category in category_map.items():
        if keyword.lower() in product:
            return category
    return "Other"

df["functional_category"] = df["product"].apply(classify)

# -----------------------------
# 4. Build summary table
# -----------------------------
# Number of genes per category
gene_counts = df.groupby("functional_category")["locus_tag"].count().reset_index()
gene_counts.columns = ["Functional category", "Number of genes"]

# Number of strains per category
strain_counts = df.groupby("functional_category")["strain_name"].nunique().reset_index()
strain_counts.columns = ["Functional category", "Number of strains"]

# Collect strain names per category
strain_list = (df.groupby("functional_category")["strain_name"]
                 .apply(lambda x: ", ".join(sorted(set(x.dropna()))))
                 .reset_index())
strain_list.columns = ["Functional category", "Strains"]

# Merge all together
summary = gene_counts.merge(strain_counts, on="Functional category")
summary = summary.merge(strain_list, on="Functional category")

# -----------------------------
# 5. Save and display
# -----------------------------
summary.to_csv("functional_category_with_strains.csv", index=False)
print(summary)

import matplotlib.pyplot as plt

# Data
traits = ['Redox / energy', 'Transport', 'Stress response', 'Other']
num_genes = [481, 190, 137, 45]
percentages = [56.4, 22.3, 16.1, 5.3]

# Create bar plot
fig, ax = plt.subplots(figsize=(8,6))
bars = ax.bar(traits, percentages, color=['#d62728', '#1f77b4', '#ff7f0e', '#2ca02c', ])

# Add data labels on top of bars
for bar, pct in zip(bars, percentages):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 1, f'{pct}%', ha='center', va='bottom', fontsize=12)

# Set labels and title
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_xlabel('Functional Trait', fontsize=12)
ax.set_title('Functional Trait Distribution of Genes', fontsize=14)
ax.set_ylim(0, 65)  # To make space for labels

plt.tight_layout()
plt.show()
