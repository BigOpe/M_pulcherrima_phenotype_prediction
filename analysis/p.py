import pandas as pd

# -----------------------------
# 1. Load contig list
# -----------------------------
contigs = []

with open(r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\contigs.txt") as f:
    for line in f:
        acc, length = line.strip().split()
        contigs.append([acc, int(length)])

df = pd.DataFrame(contigs, columns=["contig", "length"])

# -----------------------------
# 2. Extract strain prefix
# -----------------------------
# Example: CP034462.1 → CP0344
#          JAQFWG010000001.1 → JAQFWG
df["strain_prefix"] = df["contig"].str.extract(r"^([A-Z]+)")

# -----------------------------
# 3. Calculate assembly statistics
# -----------------------------
summary = (
    df.groupby("strain_prefix")
      .agg(
          Number_of_Contigs=("contig", "count"),
          Total_Length_bp=("length", "sum"),
          Longest_Contig_bp=("length", "max")
      )
      .reset_index()
      .sort_values("Total_Length_bp", ascending=False)
)

# -----------------------------
# 4. Save results
# -----------------------------
summary.to_csv("GEN_strain_contig_summary.csv", index=False)

print("Strain contig summary saved to GEN_strain_contig_summary.csv")
print(summary)



from Bio import SeqIO
import os
from collections import defaultdict

# Path to your FASTA file
fasta_file = "C:/Users/Giorgia/dna/dna2vec-legacy/gen.fasta"

# Dictionary to collect stats
stats = defaultdict(lambda: {"total_length": 0, "contigs": 0, "longest": 0})

for record in SeqIO.parse(fasta_file, "fasta"):
    seq_len = len(record.seq)
    # Use the accession prefix before the dot as strain ID
    strain_id = record.id.split(".")[0]
    
    stats[strain_id]["total_length"] += seq_len
    stats[strain_id]["contigs"] += 1
    stats[strain_id]["longest"] = max(stats[strain_id]["longest"], seq_len)

# Print summary per strain
for strain, s in stats.items():
    print(f"{strain}: {s['contigs']} contigs, total length {s['total_length']}, longest contig {s['longest']}")

from Bio import SeqIO
import os

# Input FASTA file
fasta_file = r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\gen.fasta"

# Output directory
output_dir = "C:/Users/Giorgia/dna/dna2vec-legacy/representative_genomes/"
os.makedirs(output_dir, exist_ok=True)

# Collect stats
records = list(SeqIO.parse(fasta_file, "fasta"))
total_length = sum(len(r.seq) for r in records)
contigs = len(records)
longest = max(len(r.seq) for r in records)

print(f"Assembly stats: {contigs} contigs, total length {total_length}, longest contig {longest}")

# Choose this assembly as representative (since only one file)
strain = "GEN"   # you can replace with actual strain name
output_file = os.path.join(output_dir, f"{strain}_representative.fasta")

# Rewrite headers
with open(output_file, "w") as out:
    for i, record in enumerate(records, start=1):
        record.id = f"{strain}_contig{i}"
        record.description = ""
        SeqIO.write(record, out, "fasta")

print(f"Saved representative FASTA for {strain} → {output_file}")

import pandas as pd

# Load Prokka annotation table
df = pd.read_csv(
    r"C:\Users\Giorgia\dna\dna2vec-legacy\representative_genomes\GEN_annotation\GEN.tsv",
    sep="\t"
)


# Define marker keywords (edit to match your phenotype)
marker_keywords = [
    "stress",
    "heat shock",
    "oxidative",
    "chaperone",
    "dehydrogenase",
    "transporter",
    "membrane",
    "efflux",
    "oxidoreductase"
]

# Find marker genes
marker_df = df[df["product"].str.lower().str.contains("|".join(marker_keywords), na=False)]

# Add strain name
marker_df["strain"] = "GEN"

# Save results
marker_df.to_csv("GEN_marker_genes.tsv", sep="\t", index=False)

print(marker_df.head())
print(f"Total marker genes found: {len(marker_df)}")

import pandas as pd

df = marker_df.copy()

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

df["functional_category"] = df["product"].apply(classify)

category_counts = df["functional_category"].value_counts()
print(category_counts)

def trait_map(product):
    product = str(product).lower()
    if "chaperone" in product or "dnaK" in product:
        return "Stress tolerance"
    elif "transporter" in product:
        return "Nutrient uptake"
    elif "dehydrogenase" in product:
        return "Fermentation / redox"
    else:
        return "General metabolism"

df["phenotypic_trait"] = df["product"].apply(trait_map)

df["phenotypic_trait"].value_counts()

import pandas as pd

# Load Prokka annotation
df = pd.read_csv(r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\GEN_marker_genes.tsv", sep="\t")

# Example: get ethanol tolerance-related genes
ethanol_genes = df[df['product'].str.contains("alcohol dehydrogenase", case=False, na=False)]
print(ethanol_genes[['locus_tag', 'product']])

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# -----------------------------
# 1. User Inputs / Configuration
# -----------------------------

# Path to Prokka output TSV
prokka_tsv = r"C:\Users\Giorgia\dna\dna2vec-legacy\representative_genomes\GEN_annotation\GEN.tsv"

# List of phenotypic traits and associated search keywords for annotation matching
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
# 2. Load Prokka Annotation
# -----------------------------
try:
    df = pd.read_csv(prokka_tsv, sep="\t")
    print(f"Loaded {len(df)} genes from {prokka_tsv}")
except FileNotFoundError:
    print(f"Error: {prokka_tsv} not found!")
    exit()

# -----------------------------
# 3. Initialize Gene–Trait Matrix
# -----------------------------
gene_trait_matrix = pd.DataFrame()
gene_trait_matrix['locus_tag'] = df['locus_tag']
gene_trait_matrix['product'] = df['product']

# Initialize trait columns with 0
for trait in traits_keywords.keys():
    gene_trait_matrix[trait] = 0

# -----------------------------
# 4. Map Genes to Traits
# -----------------------------
for trait, keywords in traits_keywords.items():
    mask = df['product'].str.contains('|'.join(keywords), case=False, na=False)
    gene_trait_matrix.loc[mask, trait] = 1

# -----------------------------
# 5. Save the Matrix
# -----------------------------
output_file = "GEN_gene_trait_matrix.csv"
gene_trait_matrix.to_csv(output_file, index=False)
print(f"Gene–trait matrix saved to {output_file}")

# -----------------------------
# 6. Optional: Trait Summary
# -----------------------------
trait_summary = gene_trait_matrix.iloc[:, 2:].sum().sort_values(ascending=False)
print("\nNumber of genes per trait:")
print(trait_summary)

# -----------------------------
# 7. Add functional categories
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

df["functional_category"] = df["product"].apply(classify)

# -----------------------------
# 8. Melt gene–trait matrix to long format
# -----------------------------
trait_df = gene_trait_matrix.melt(
    id_vars=['locus_tag', 'product'],
    var_name='trait',
    value_name='has_trait'
)
trait_df = trait_df[trait_df['has_trait'] == 1].drop(columns='has_trait')

# Merge annotation with traits
merged_df = pd.merge(df, trait_df, on="locus_tag", how="inner")

# -----------------------------
# 9. Build a contingency table: traits vs functional categories
# -----------------------------
all_traits = list(traits_keywords.keys())
heatmap_df = merged_df.groupby(['functional_category', 'trait']).size().unstack(fill_value=0)

# Ensure all traits are included
for trait in all_traits:
    if trait not in heatmap_df.columns:
        heatmap_df[trait] = 0
heatmap_df = heatmap_df[all_traits]  # reorder columns

# -----------------------------
# 10. Plot heatmap
# -----------------------------
plt.figure(figsize=(12,6))
sns.heatmap(heatmap_df, annot=True, fmt="d", cmap="YlGnBu")
plt.title("Gene counts per Trait vs Functional Category")
plt.ylabel("Functional Category")
plt.xlabel("Trait")
plt.tight_layout()
plt.show()

# -----------------------------
# 11. Save heatmap data
# -----------------------------
heatmap_df.to_csv("GEN_trait_category_counts.csv")
print("Saved trait-category counts to GEN_trait_category_counts.csv")
