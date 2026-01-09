import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. User Inputs / Configuration
# -----------------------------
prokka_tsv = r"C:\Users\Giorgia\dna\dna2vec-legacy\representative_genomes\GEN_annotation\GEN.tsv"

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
csv_file = r"C:\Users\Giorgia\dna\GEN_analysis_results\GEN_trait_category_percentage.csv"
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
heatmap_file = r"C:\Users\Giorgia\dna\GEN_analysis_results\GEN_trait_category_percentage.csv"
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

