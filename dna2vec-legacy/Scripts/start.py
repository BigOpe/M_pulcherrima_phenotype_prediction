import os
import pandas as pd

# -----------------------------
# 1. Configuration
# -----------------------------
prokka_parent_dir = r"C:\Users\Giorgia\dna"  # Parent folder containing all prokka_* folders
output_file = r"C:\Users\Giorgia\dna\all_strains_full_gene_table.csv"

# Functional classification function
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

# -----------------------------
# 2. Collect Prokka TSV files
# -----------------------------
all_gene_dfs = []

for folder_name in os.listdir(prokka_parent_dir):
    if folder_name.startswith("prokka_"):
        folder_path = os.path.join(prokka_parent_dir, folder_name)
        if os.path.isdir(folder_path):
            # Find the TSV file in the folder
            tsv_files = [f for f in os.listdir(folder_path) if f.endswith(".tsv")]
            if tsv_files:
                tsv_file = tsv_files[0]  # Usually only one per folder
                tsv_path = os.path.join(folder_path, tsv_file)
                
                # Load TSV
                df = pd.read_csv(tsv_path, sep="\t")
                
                # Add strain info from folder name
                strain_name = folder_name.replace("prokka_", "")
                df["strain"] = strain_name
                
                all_gene_dfs.append(df)
                print(f"Loaded {len(df)} genes from {strain_name}")
            else:
                print(f"No TSV found in {folder_name}")

# -----------------------------
# 3. Concatenate all TSVs
# -----------------------------
if all_gene_dfs:
    full_gene_table = pd.concat(all_gene_dfs, ignore_index=True)
    print(f"\nTotal genes collected: {len(full_gene_table)}")
else:
    raise ValueError("No Prokka TSV files found!")

# -----------------------------
# 4. Classify functional categories
# -----------------------------
full_gene_table["functional_category"] = full_gene_table["product"].apply(classify)

# -----------------------------
# 5. Save full table
# -----------------------------
full_gene_table.to_csv(output_file, index=False)
print(f"Saved combined table to {output_file}")

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. User Configuration
# -----------------------------
# Folder containing all Prokka annotation TSVs
prokka_folder = r"C:\Users\Giorgia\dna"

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
# 2. Load all TSVs and add strain column
# -----------------------------
all_genes = []

for folder in os.listdir(prokka_folder):
    folder_path = os.path.join(prokka_folder, folder)
    if os.path.isdir(folder_path) and folder.startswith("prokka_"):
        strain_id = folder.replace("prokka_", "")
        tsv_files = [f for f in os.listdir(folder_path) if f.endswith(".tsv")]
        for tsv in tsv_files:
            tsv_path = os.path.join(folder_path, tsv)
            df = pd.read_csv(tsv_path, sep="\t")
            df['strain'] = strain_id
            all_genes.append(df)

# Combine into a single DataFrame
if not all_genes:
    raise ValueError("No Prokka TSV files found!")
full_gene_table = pd.concat(all_genes, ignore_index=True)
print(f"Total genes collected: {len(full_gene_table)}")

# -----------------------------
# 3. Map genes to traits
# -----------------------------
for trait, keywords in traits_keywords.items():
    mask = full_gene_table['product'].str.contains('|'.join(keywords), case=False, na=False)
    full_gene_table[trait] = 0
    full_gene_table.loc[mask, trait] = 1

# -----------------------------
# 4. Add functional categories
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

full_gene_table["functional_category"] = full_gene_table["product"].apply(classify)

# -----------------------------
# 5. Save gene–trait matrix
# -----------------------------
output_matrix_file = os.path.join(prokka_folder, "all_strains_gene_trait_matrix.csv")
cols_to_save = ['locus_tag', 'product', 'strain'] + list(traits_keywords.keys())
full_gene_table[cols_to_save].to_csv(output_matrix_file, index=False)
print(f"Gene–trait matrix saved to {output_matrix_file}")

# -----------------------------
# 6. Build trait vs functional category heatmap
# -----------------------------
# Melt to long format
trait_df = full_gene_table.melt(
    id_vars=['locus_tag', 'product', 'functional_category'],
    var_name='trait',
    value_name='has_trait'
)
trait_df = trait_df[trait_df['has_trait'] == 1].drop(columns='has_trait')

# Group by functional category and trait
heatmap_df = trait_df.groupby(['functional_category', 'trait']).size().unstack(fill_value=0)

# Ensure all traits are present
for trait in traits_keywords.keys():
    if trait not in heatmap_df.columns:
        heatmap_df[trait] = 0
heatmap_df = heatmap_df[list(traits_keywords.keys())]  # reorder

# Plot heatmap
plt.figure(figsize=(12,6))
sns.heatmap(heatmap_df, annot=True, fmt="d", cmap="YlGnBu")
plt.title("Gene counts per Trait vs Functional Category")
plt.ylabel("Functional Category")
plt.xlabel("Trait")
plt.tight_layout()
plt.show()

# Save heatmap data
heatmap_output_file = os.path.join(prokka_folder, "all_strains_trait_category_counts.csv")
heatmap_df.to_csv(heatmap_output_file)
print(f"Saved trait-category counts to {heatmap_output_file}")

import os
import pandas as pd

# -----------------------------
# 1. User Configuration
# -----------------------------
# Folder containing all Prokka annotation TSVs
prokka_folder = r"C:\Users\Giorgia\dna"

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
