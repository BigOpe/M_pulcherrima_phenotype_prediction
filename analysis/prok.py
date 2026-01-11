import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. User Inputs / Configuration
# -----------------------------
prokka_tsv = r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\GEN.tsv"

# Define your traits and keywords for annotation matching
traits_keywords = {
    "Ethanol_Tolerance": ["alcohol dehydrogenase", "adh", "stress response"],
    "Copper_Resistance": ["cop", "metallothionein", "copper"],
    "Sugar_Utilization": ["glycosidase", "PTS transporter", "sugar transport"],
    "H2S": ["cysteine desulfhydrase", "sulfite reductase"],
    "Protease": ["serine protease", "metalloprotease"],
    "Lipase": ["lipase"],
    "Growth_Rate": ["ribosomal protein", "ATP synthase", "energy metabolism"],
    "Glucosidase": ["glucosidase"],
    "SO2_Resistance": ["sulfite reductase", "stress response"],
    "Esterase": ["esterase"],
    "Gluconic_Acid_Utilization": ["gluconate kinase", "gluconate transporter"]
}

# -----------------------------
# 2. Load Prokka Annotation
# -----------------------------
try:
    df = pd.read_csv(prokka_tsv, sep="\t")
    print(f"Loaded {len(df)} genes from {prokka_tsv}")
except FileNotFoundError:
    print(f"Error: {prokka_tsv} not found!")
    exit()

# Ensure 'product' column exists
if 'product' not in df.columns:
    df['product'] = ""  # fill with empty strings if missing

# -----------------------------
# 3. Map Genes to Traits
# -----------------------------
gene_trait_matrix = df[['locus_tag', 'product']].copy()
for trait, keywords in traits_keywords.items():
    mask = df['product'].str.contains('|'.join(keywords), case=False, na=False)
    gene_trait_matrix[trait] = mask.astype(int)

# Save gene–trait matrix
gene_trait_matrix.to_csv("GEN_gene_trait_matrix.csv", index=False)
print("Gene–trait matrix saved to GEN_gene_trait_matrix.csv")

# -----------------------------
# 4. Add Functional Categories
# -----------------------------
def classify(product):
    product = str(product).lower()
    if "transporter" in product or "permease" in product:
        return "Transport"
    elif "chaperone" in product or "heat shock" in product:
        return "Stress response"
    elif "dehydrogenase" in product or "oxidoreductase" in product:
        return "Redox / energy"
    elif "kinase" in product or "regulator" in product:
        return "Regulation"
    elif "synthase" in product or "ligase" in product:
        return "Metabolism"
    else:
        return "Other"

df['functional_category'] = df['product'].apply(classify)

# -----------------------------
# 5. Merge gene annotation with traits
# -----------------------------
# Convert wide gene_trait_matrix to long format for merging
trait_long = gene_trait_matrix.melt(
    id_vars=['locus_tag', 'product'],
    var_name='trait',
    value_name='presence'
)
trait_long = trait_long[trait_long['presence'] == 1].drop(columns='presence')

# Merge with functional categories
merged_df = pd.merge(df, trait_long, on=['locus_tag', 'product'], how='inner')

# -----------------------------
# 6. Build Heatmap: Functional Category vs Trait
# -----------------------------
heatmap_df = merged_df.groupby(['functional_category', 'trait']).size().unstack(fill_value=0)

plt.figure(figsize=(12, 6))
sns.heatmap(heatmap_df, annot=True, fmt='d', cmap='YlGnBu')
plt.title("Gene counts per Trait vs Functional Category")
plt.ylabel("Functional Category")
plt.xlabel("Trait")
plt.tight_layout()
plt.show()

# Save heatmap data
heatmap_df.to_csv("GEN_trait_category_counts.csv")
print("Saved trait-category counts to GEN_trait_category_counts.csv")

# -----------------------------
# 7. Save top genes per trait
# -----------------------------
for trait_name in merged_df['trait'].unique():
    trait_genes = merged_df.loc[
        merged_df['trait'] == trait_name,
        ['locus_tag', 'product', 'functional_category']
    ]
    output_csv = f"top_genes_{trait_name}.csv"
    trait_genes.to_csv(output_csv, index=False)
    print(f"Saved {len(trait_genes)} genes for trait '{trait_name}' to {output_csv}")

# -----------------------------
# 8. Trait Summary
# -----------------------------
trait_counts = merged_df['trait'].value_counts()
print("\nGene counts per trait:")
print(trait_counts)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------------------
# 1. Load merged annotation + gene-trait data
# -------------------------------
merged_file = "GEN_gene_trait_matrix.csv"  # Your existing gene-trait matrix
df_trait = pd.read_csv(merged_file)

# -------------------------------
# 2. Ensure functional_category exists
# If not, add it (optional)
# -------------------------------
if "functional_category" not in df_trait.columns:
    def classify(product):
        product = str(product).lower()
        if "transporter" in product or "permease" in product:
            return "Transport"
        elif "chaperone" in product or "heat shock" in product:
            return "Stress response"
        elif "dehydrogenase" in product or "oxidoreductase" in product:
            return "Redox / energy"
        elif "kinase" in product or "regulator" in product:
            return "Regulation"
        elif "synthase" in product or "ligase" in product:
            return "Metabolism"
        else:
            return "Other"
    df_trait["functional_category"] = df_trait["product"].apply(classify)

# -------------------------------
# 3. Melt the gene-trait matrix
# -------------------------------
traits = df_trait.columns[2:]  # Assumes first 2 columns are locus_tag and product
melted = df_trait.melt(id_vars=["locus_tag", "product", "functional_category"],
                       value_vars=traits,
                       var_name="trait",
                       value_name="presence")
# Keep only genes that are associated with a trait
melted = melted[melted['presence'] == 1]

# -------------------------------
# 4. Contingency table: functional category vs trait
# -------------------------------
heatmap_df = pd.crosstab(melted['functional_category'], melted['trait'])

# Convert to percentages per trait
heatmap_pct = heatmap_df.div(heatmap_df.sum(axis=0), axis=1) * 100

# -------------------------------
# 5. Plot heatmap
# -------------------------------
plt.figure(figsize=(12,6))
sns.heatmap(heatmap_pct, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': '% of genes'})
plt.title("Percentage of Genes per Trait vs Functional Category")
plt.ylabel("Functional Category")
plt.xlabel("Trait")
plt.tight_layout()
plt.show()

# -------------------------------
# 6. Save percentage table
# -------------------------------
heatmap_pct.to_csv("GEN_trait_category_percentage.csv")
print("Saved trait-category percentages to GEN_trait_category_percentage.csv")

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------
# 1. Load the CSV
# -----------------------------
csv_file = r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\GEN_trait_category_percentage.csv"
df = pd.read_csv(csv_file, index_col=0)  # Functional categories as index
print("Loaded data:")
print(df.head())

# -----------------------------
# 2. Compute total counts per trait
# -----------------------------
total_per_trait = df.sum(axis=0)
print("\nTotal genes per trait (sum of all categories):")
print(total_per_trait)

# -----------------------------
# 3. Compute fraction (%) of each functional category per trait
# -----------------------------
fraction_df = df.div(df.sum(axis=0), axis=1) * 100  # Convert to percentages
print("\nFraction (%) of genes per functional category for each trait:")
print(fraction_df)

# -----------------------------
# 4. Determine top 3 functional categories per trait
# -----------------------------
top_categories_list = []

for trait in fraction_df.columns:
    sorted_categories = fraction_df[trait].sort_values(ascending=False)
    top3 = sorted_categories.head(3)
    cumulative = 0
    for rank, (category, pct) in enumerate(top3.items(), start=1):
        cumulative += pct
        top_categories_list.append({
            "Trait": trait,
            "Rank": rank,
            "Functional_Category": category,
            "Percentage": pct,
            "Cumulative_Percentage": cumulative
        })

top3_df = pd.DataFrame(top_categories_list)
print("\nTop 3 functional categories per trait with cumulative percentage:")
print(top3_df)

# -----------------------------
# 5. Save top 3 summary table
# -----------------------------
summary_file = "GEN_trait_top3_categories_cumulative.csv"
top3_df.to_csv(summary_file, index=False)
print(f"\nSaved top 3 functional categories per trait (with cumulative %) to {summary_file}")

# -----------------------------
# 6. Plot stacked bar chart
# -----------------------------
plt.figure(figsize=(12,6))
fraction_df.T.plot(kind='bar', stacked=True, colormap='tab20', figsize=(12,6))
plt.ylabel("Percentage of genes")
plt.xlabel("Trait")
plt.title("Functional Category Distribution per Trait")
plt.legend(title="Functional Category", bbox_to_anchor=(1.05, 1))
plt.tight_layout()
plt.show()

import pandas as pd

# -------------------------------
# 1. Load trait-category percentage data
# -------------------------------
heatmap_file = r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\GEN_trait_category_percentage.csv"
heatmap_df = pd.read_csv(heatmap_file, index_col=0)  # rows = functional categories

# -------------------------------
# 2. Compute top 3 functional categories per trait
# -------------------------------
top3_list = []

for trait in heatmap_df.columns:
    # Sort functional categories for this trait by percentage descending
    sorted_df = heatmap_df[[trait]].sort_values(by=trait, ascending=False)
    sorted_df = sorted_df.reset_index()  # reset index to make 'functional_category' a column
    sorted_df['Trait'] = trait
    sorted_df['Rank'] = range(1, len(sorted_df)+1)  # rank each row
    top3_list.append(sorted_df.head(3))  # keep only top 3

# Combine all traits into a single DataFrame
top3_categories = pd.concat(top3_list, ignore_index=True)

# Rename columns for clarity
top3_categories.rename(columns={'functional_category': 'Functional_Category',
                                trait: 'Percentage'}, inplace=True)

# Compute cumulative percentage per trait
top3_categories['Cumulative_Percentage'] = top3_categories.groupby('Trait')['Percentage'].cumsum()

# -------------------------------
# 3. Save to CSV
# -------------------------------
output_file = "GEN_trait_top3_categories_cumulative.csv"
top3_categories.to_csv(output_file, index=False)
print(f"Saved top 3 functional categories per trait (with cumulative %) to {output_file}")

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


import pandas as pd

# Load combined gene–trait / Prokka annotation table
df = pd.read_csv(r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\all_strains_full_gene_table.csv")


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
print("Columns in df:", df.columns.tolist())
print(df.head())

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
# =====================================
# 0. Imports
# =====================================
import os
import pandas as pd
import matplotlib.pyplot as plt
import pingouin as pg


# =====================================
# 1. Load gene–trait matrix
# =====================================
gene_df = pd.read_csv(
    r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\all_strains_full_gene_table.csv")

print("Gene DF columns:", gene_df.columns.tolist())


# =====================================
# 2. Prefix → Strain mapping
# =====================================
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

# Extract prefix (letters at start of strain/contig)
gene_df["prefix"] = gene_df["strain"].astype(str).str.extract(r"^([A-Z]+)")
gene_df["strain_name"] = gene_df["prefix"].map(prefix_map)


# =====================================
# 3. Functional category classification
# =====================================
category_map = {
    "dehydrogenase": "Redox / Energy metabolism",
    "oxidoreductase": "Redox / Energy metabolism",
    "reductase": "Redox / Energy metabolism",
    "alcohol": "Ethanol metabolism",
    "transport": "Transporters",
    "permease": "Transporters",
    "abc": "Transporters",
    "cell wall": "Cell wall organization",
    "mannoprotein": "Cell wall organization",
    "glucan": "Cell wall organization",
    "stress": "Stress response",
    "heat shock": "Stress response",
    "hsp": "Stress response",
    "lipase": "Lipid metabolism",
    "esterase": "Lipid metabolism",
    "protease": "Protein degradation",
    "ubiquitin": "Protein turnover",
    "dna": "Indeterminate / housekeeping",
    "rna": "Indeterminate / housekeeping"
}

def classify_gene(product):
    product = str(product).lower()
    for keyword, category in category_map.items():
        if keyword in product:
            return category
    return "Other"

gene_df["functional_category"] = gene_df["product"].apply(classify_gene)


# =====================================
# 4. Summary tables
# =====================================
gene_counts = (
    gene_df.groupby("functional_category")["locus_tag"]
    .count()
    .reset_index(name="Number of genes")
)

strain_counts = (
    gene_df.groupby("functional_category")["strain_name"]
    .nunique()
    .reset_index(name="Number of strains")
)

strain_list = (
    gene_df.groupby("functional_category")["strain_name"]
    .apply(lambda x: ", ".join(sorted(set(x.dropna()))))
    .reset_index(name="Strains")
)

summary_df = (
    gene_counts
    .merge(strain_counts, on="functional_category")
    .merge(strain_list, on="functional_category")
)

summary_df.rename(columns={"functional_category": "Functional category"}, inplace=True)

summary_df.to_csv("functional_category_with_strains.csv", index=False)
print(summary_df)


# =====================================
# 5. Plot functional distribution
# =====================================
plot_df = summary_df.sort_values("Number of genes", ascending=False)

plt.figure(figsize=(8, 6))
plt.bar(plot_df["Functional category"], plot_df["Number of genes"])
plt.xticks(rotation=45, ha="right")
plt.ylabel("Number of genes")
plt.title("Functional Trait Distribution of Genes")
plt.tight_layout()
plt.show()


# =====================================
# 6. Phenotype analysis (Cronbach alpha)
# =====================================
pheno_df = pd.read_csv(
    r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\phen.csv"
)

pheno_df = pheno_df.set_index("Strain")
pheno_numeric = pheno_df.drop(columns=["Notes"], errors="ignore")

pheno_numeric = pheno_numeric.apply(pd.to_numeric, errors="coerce")
pheno_numeric = pheno_numeric.fillna(pheno_numeric.mean())

alpha, ci = pg.cronbach_alpha(pheno_numeric)

print(f"Cronbach’s alpha: {alpha:.3f}")
print(f"95% CI: {ci}")



import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. User Configuration
# -----------------------------
# Folder containing all Prokka annotation TSVs
prokka_folder = r"C:\Users\Giorgia\Desktop\phd-ope\dna"

# List of phenotypic traits and associated keywords
traits_keywords = {
    "Ethanol_Tolerance": ["alcohol dehydrogenase", "adh", "stress response"],
    "Killer_Activity": ["killer toxin", "k1", "k2", "secreted protein"],
    "Pigment_Intensity": ["melanin", "carotenoid"],
    "Copper_Resistance": ["cop", "metallothionein", "copper"],
    "Sugar_Utilization": ["glycosidase", "PTS transporter", "sugar transport"],
    "H2S": ["cysteine desulfhydrase", "sulfite reductase"],
    "Protease": ["serine protease", "metalloprotease"],
    "Lipase": ["lipase"],
    "Growth_Rate": ["ribosomal protein", "ATP synthase", "energy metabolism"],
    "Glucosidase": ["glucosidase"],
    "SO2_Resistance": ["sulfite reductase", "stress response"],
    "Esterase": ["esterase"],
    "Gluconic_Acid_Utilization": ["gluconate kinase", "gluconate transporter"]
}


# -----------------------------
# 4. Add functional categories
# -----------------------------

# Group by functional category and trait

import os
import pandas as pd

# -----------------------------
# 1. User Configuration
# -----------------------------
# Folder containing all Prokka annotation TSVs
prokka_folder = r"C:\Users\Giorgia\Desktop\phd-ope"

# List of phenotypic traits and associated keywords
traits_keywords = {
    "Ethanol_Tolerance": ["alcohol dehydrogenase", "adh", "stress response"],
    "Killer_Activity": ["killer toxin", "k1", "k2", "secreted protein"],
    "Pigment_Intensity": ["melanin", "carotenoid"],
    "Copper_Resistance": ["cop", "metallothionein", "copper"],
    "Sugar_Utilization": ["glycosidase", "PTS transporter", "sugar transport"],
    "H2S": ["cysteine desulfhydrase", "sulfite reductase"],
    "Protease": ["serine protease", "metalloprotease"],
    "Lipase": ["lipase"],
    "Growth_Rate": ["ribosomal protein", "ATP synthase", "energy metabolism"],
    "Glucosidase": ["glucosidase"],
    "SO2_Resistance": ["sulfite reductase", "stress response"],
    "Esterase": ["esterase"],
    "Gluconic_Acid_Utilization": ["gluconate kinase", "gluconate transporter"]
}

# -----------------------------
# 2. Process each strain individually
# -----------------------------
for folder in os.listdir(prokka_folder):
    folder_path = os.path.join(prokka_folder, folder)
    if os.path.isdir(folder_path) and folder.startswith("prokka_"):
        strain_id = folder.replace("prokka_", "")
        tsv_files = [f for f in os.listdir(folder_path) if f.endswith(".tsv")]
        all_genes = []
        for tsv in tsv_files:
            tsv_path = os.path.join(folder_path, tsv)
            df = pd.read_csv(tsv_path, sep="\t")
            df['strain'] = strain_id
            all_genes.append(df)
        if not all_genes:
            continue
        full_gene_table = pd.concat(all_genes, ignore_index=True)
        
        # Map genes to traits
        for trait, keywords in traits_keywords.items():
            mask = full_gene_table['product'].str.contains('|'.join(keywords), case=False, na=False)
            full_gene_table[trait] = 0
            full_gene_table.loc[mask, trait] = 1
        
        # Add functional categories
        def classify(product):
            product = str(product).lower()
            if "transporter" in product or "permease" in product:
                return "Transport"
            elif "chaperone" in product or "heat shock" in product:
                return "Stress response"
            elif "dehydrogenase" in product or "oxidoreductase" in product:
                return "Redox / energy"
            elif "kinase" in product or "regulator" in product:
                return "Regulation"
            elif "synthase" in product or "ligase" in product:
                return "Metabolism"
            else:
                return "Other"
        
        full_gene_table["functional_category"] = full_gene_table["product"].apply(classify)
        
        # Save per-strain gene–trait matrix
        output_file = os.path.join(prokka_folder, f"{strain_id}_gene_trait_matrix.csv")
        cols_to_save = ['locus_tag', 'product', 'strain'] + list(traits_keywords.keys())
        full_gene_table[cols_to_save].to_csv(output_file, index=False)
        print(f"Saved gene–trait matrix for strain {strain_id} → {output_file}")
        
