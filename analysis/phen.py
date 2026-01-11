import pandas as pd

# -----------------------------
# 1. Load contig list
# -----------------------------
contigs = []

with open(r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\contigs.txt") as f:
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

# Path to your FASTA file
fasta_file = "C:/Users/Giorgia/dna/dna2vec-legacy/gen.fasta"

# Parse the FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    print(record.id, len(record.seq))

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
fasta_file = "C:/Users/Giorgia/dna/dna2vec-legacy/gen.fasta"

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
df = pd.read_csv(r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\GEN.tsv", sep="\t")

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

# Aggregating contig lengths by strain prefix and summarizing contig statistics

import re
import pandas as pd
from collections import defaultdict

# Raw input string (truncated for brevity in this example)
raw_data = """CP034462.1 687369 CP034461.1 1422419 CP034460.1 1923762 JAQFWG010000001.1 1304897 JAQFWG010000002.1 896825 JAQFWG010000003.1 120734 JAJMIJ010000354.1 35377 JAJMIJ010000346.1 39091 JAJMIJ010000347.1 38308 JAJMIJ010000348.1 36485 NISE01000353.1 7362 NISE01000345.1 7959 NISE01000346.1 7930 NISE01000347.1 7793 JAKTYT010000001.1 46026 JAKTYT010000002.1 42237 JAKTYT010000003.1 37844 JAKTYT010000005.1 32328 QBLL01000001.1 106026 QBLL01000002.1 32332 QBLL01000003.1 31415 QBLL01000004.1 30875 JAJMIJ010000349.1 36371 JAJMIJ010000350.1 36295 JAJMIJ010000351.1 36279 JAJMIJ010000352.1 35701 JAJMIJ010000353.1 35408 JAJMIJ010000345.1 39768 JAJMIJ010000355.1 34868 JAJMIJ010000356.1 34138 JAJMIJ010000357.1 33098 JAJMIJ010000358.1 32868 JAJMIJ010000360.1 30968 JAJMIJ010000361.1 30850 JAJMIJ010000362.1 30681 JAJMIJ010000344.1 40345 JAJMIJ010000343.1 40587"""

# Step 1: Parse contig ID and length pairs
tokens = raw_data.strip().split()
contig_data = [(tokens[i], int(tokens[i+1])) for i in range(0, len(tokens), 2)]

# Step 2: Group by strain prefix
strain_stats = defaultdict(list)

for contig_id, length in contig_data:
    # Extract strain prefix (everything before the first digit sequence)
    match = re.match(r"([A-Z]+[A-Z0-9]*)(\d{3,})", contig_id.replace('.', ''))
    if match:
        strain_prefix = match.group(1)
    else:
        strain_prefix = contig_id.split('.')[0]
    strain_stats[strain_prefix].append(length)

# Step 3: Summarize per strain
summary = []
for strain, lengths in strain_stats.items():
    summary.append({
        "Strain": strain,
        "Number of Contigs": len(lengths),
        "Total Length (bp)": sum(lengths),
        "Longest Contig (bp)": max(lengths)
    })

# Step 4: Create DataFrame and sort
df_summary = pd.DataFrame(summary)
df_summary = df_summary.sort_values(by="Total Length (bp)", ascending=False)

# Save to file
output_path = "/mnt/data/strain_contig_summary.csv"
df_summary.to_csv(output_path, index=False)

print("Aggregated contig lengths by strain and saved summary to strain_contig_summary.csv")



import pandas as pd

# Step 1: Create contig-to-strain mapping
contig_to_strain = {
    # CP
    "CP034462.1": "CP",
    "CP034461.1": "CP",
    "CP034460.1": "CP",
    "CP034459.1": "CP",
    "CP034458.1": "CP",
    "CP034457.1": "CP",
    "CP034456.1": "CP",
    
    # JAQFWG
    "JAQFWG010000001.1": "JAQFWG",
    "JAQFWG010000002.1": "JAQFWG",
    "JAQFWG010000003.1": "JAQFWG",
    "JAQFWG010000006.1": "JAQFWG",
    "JAQFWG010000007.1": "JAQFWG",
    "JAQFWG010000008.1": "JAQFWG",
    "JAQFWG010000009.1": "JAQFWG",
    "JAQFWG010000011.1": "JAQFWG",
    "JAQFWG010000012.1": "JAQFWG",
    "JAQFWG010000013.1": "JAQFWG",
    "JAQFWG010000014.1": "JAQFWG",
    "JAQFWG010000015.1": "JAQFWG",
    "JAQFWG010000016.1": "JAQFWG",
    "JAQFWG010000017.1": "JAQFWG",
    "JAQFWG010000018.1": "JAQFWG",
    "JAQFWG010000019.1": "JAQFWG",
    "JAQFWG010000020.1": "JAQFWG",
    "JAQFWG010000021.1": "JAQFWG",
    "JAQFWG010000022.1": "JAQFWG",
    "JAQFWG010000023.1": "JAQFWG",
    "JAQFWG010000024.1": "JAQFWG",
    
    # JAJMIJ
    "JAJMIJ010000343.1": "JAJMIJ",
    "JAJMIJ010000344.1": "JAJMIJ",
    "JAJMIJ010000345.1": "JAJMIJ",
    "JAJMIJ010000346.1": "JAJMIJ",
    "JAJMIJ010000347.1": "JAJMIJ",
    "JAJMIJ010000348.1": "JAJMIJ",
    "JAJMIJ010000349.1": "JAJMIJ",
    "JAJMIJ010000350.1": "JAJMIJ",
    "JAJMIJ010000351.1": "JAJMIJ",
    "JAJMIJ010000352.1": "JAJMIJ",
    "JAJMIJ010000353.1": "JAJMIJ",
    "JAJMIJ010000354.1": "JAJMIJ",
    "JAJMIJ010000355.1": "JAJMIJ",
    "JAJMIJ010000356.1": "JAJMIJ",
    "JAJMIJ010000357.1": "JAJMIJ",
    "JAJMIJ010000358.1": "JAJMIJ",
    "JAJMIJ010000360.1": "JAJMIJ",
    "JAJMIJ010000361.1": "JAJMIJ",
    "JAJMIJ010000362.1": "JAJMIJ",
    
    # NISE
    "NISE01000342.1": "NISE",
    "NISE01000343.1": "NISE",
    "NISE01000344.1": "NISE",
    "NISE01000345.1": "NISE",
    "NISE01000346.1": "NISE",
    "NISE01000347.1": "NISE",
    "NISE01000348.1": "NISE",
    "NISE01000349.1": "NISE",
    "NISE01000350.1": "NISE",
    "NISE01000351.1": "NISE",
    "NISE01000352.1": "NISE",
    "NISE01000353.1": "NISE",
    "NISE01000354.1": "NISE",
    "NISE01000355.1": "NISE",
    "NISE01000356.1": "NISE",
    "NISE01000357.1": "NISE",
    "NISE01000358.1": "NISE",
    "NISE01000359.1": "NISE",
    "NISE01000360.1": "NISE",
    "NISE01000361.1": "NISE",
    
    # JAKTYT
    "JAKTYT010000001.1": "JAKTYT",
    "JAKTYT010000002.1": "JAKTYT",
    "JAKTYT010000003.1": "JAKTYT",
    "JAKTYT010000005.1": "JAKTYT",
    "JAKTYT010000007.1": "JAKTYT",
    "JAKTYT010000008.1": "JAKTYT",
    "JAKTYT010000009.1": "JAKTYT",
    "JAKTYT010000010.1": "JAKTYT",
    "JAKTYT010000011.1": "JAKTYT",
    "JAKTYT010000012.1": "JAKTYT",
    "JAKTYT010000013.1": "JAKTYT",
    "JAKTYT010000015.1": "JAKTYT",
    "JAKTYT010000017.1": "JAKTYT",
    "JAKTYT010000018.1": "JAKTYT",
    "JAKTYT010000019.1": "JAKTYT",
    "JAKTYT010000020.1": "JAKTYT",
    
    # QBLL
    "QBLL01000001.1": "QBLL",
    "QBLL01000002.1": "QBLL",
    "QBLL01000003.1": "QBLL",
    "QBLL01000004.1": "QBLL",
    "QBLL01000005.1": "QBLL",
    "QBLL01000006.1": "QBLL",
    "QBLL01000007.1": "QBLL",
    "QBLL01000008.1": "QBLL",
    "QBLL01000009.1": "QBLL",
    "QBLL01000010.1": "QBLL",
    "QBLL01000011.1": "QBLL",
    "QBLL01000012.1": "QBLL",
    "QBLL01000013.1": "QBLL",
    "QBLL01000014.1": "QBLL",
    "QBLL01000015.1": "QBLL",
    "QBLL01000016.1": "QBLL",
    "QBLL01000017.1": "QBLL",
    "QBLL01000018.1": "QBLL",
    "QBLL01000019.1": "QBLL",
    "QBLL01000020.1": "QBLL",
    
    # JAJMIQ
    "JAJMIQ010000001.1": "JAJMIQ",
    "JAJMIQ010000002.1": "JAJMIQ",
    "JAJMIQ010000003.1": "JAJMIQ",
    "JAJMIQ010000004.1": "JAJMIQ",
    "JAJMIQ010000005.1": "JAJMIQ",
    "JAJMIQ010000006.1": "JAJMIQ",
    "JAJMIQ010000007.1": "JAJMIQ",
    "JAJMIQ010000008.1": "JAJMIQ",
    "JAJMIQ010000009.1": "JAJMIQ",
    "JAJMIQ010000010.1": "JAJMIQ",
    "JAJMIQ010000011.1": "JAJMIQ",
    "JAJMIQ010000012.1": "JAJMIQ",
    "JAJMIQ010000013.1": "JAJMIQ",
    "JAJMIQ010000014.1": "JAJMIQ",
    "JAJMIQ010000015.1": "JAJMIQ",
    "JAJMIQ010000016.1": "JAJMIQ",
    "JAJMIQ010000017.1": "JAJMIQ",
    "JAJMIQ010000018.1": "JAJMIQ",
    "JAJMIQ010000019.1": "JAJMIQ",
    "JAJMIQ010000020.1": "JAJMIQ",
    
    # JAJMHS
    "JAJMHS010000001.1": "JAJMHS",
    "JAJMHS010000002.1": "JAJMHS",
    "JAJMHS010000003.1": "JAJMHS",
    "JAJMHS010000004.1": "JAJMHS",
    "JAJMHS010000005.1": "JAJMHS",
    "JAJMHS010000006.1": "JAJMHS",
    "JAJMHS010000007.1": "JAJMHS",
    "JAJMHS010000008.1": "JAJMHS",
    "JAJMHS010000009.1": "JAJMHS",
    "JAJMHS010000010.1": "JAJMHS",
    "JAJMHS010000011.1": "JAJMHS",
    "JAJMHS010000012.1": "JAJMHS",
    "JAJMHS010000013.1": "JAJMHS",
    "JAJMHS010000014.1": "JAJMHS",
    "JAJMHS010000015.1": "JAJMHS",
    "JAJMHS010000016.1": "JAJMHS",
    "JAJMHS010000017.1": "JAJMHS",
    "JAJMHS010000018.1": "JAJMHS",
    "JAJMHS010000019.1": "JAJMHS",
    "JAJMHS010000020.1": "JAJMHS",
    
    # JAKTYM
    "JAKTYM010000001.1": "JAKTYM",
    "JAKTYM010000002.1": "JAKTYM",
    "JAKTYM010000003.1": "JAKTYM",
    "JAKTYM010000004.1": "JAKTYM",
    "JAKTYM010000005.1": "JAKTYM",
    "JAKTYM010000006.1": "JAKTYM",
    "JAKTYM010000007.1": "JAKTYM",
    "JAKTYM010000008.1": "JAKTYM",
    "JAKTYM010000009.1": "JAKTYM",
    "JAKTYM010000010.1": "JAKTYM",
    "JAKTYM010000011.1": "JAKTYM",
    "JAKTYM010000012.1": "JAKTYM",
    "JAKTYM010000013.1": "JAKTYM",
    "JAKTYM010000014.1": "JAKTYM",
    "JAKTYM010000015.1": "JAKTYM",
    "JAKTYM010000016.1": "JAKTYM",
    "JAKTYM010000017.1": "JAKTYM",
    "JAKTYM010000018.1": "JAKTYM",
    "JAKTYM010000019.1": "JAKTYM",
    "JAKTYM010000020.1": "JAKTYM",
    
    # JACBPP
    "JACBPP010000001.1": "JACBPP",
    "JACBPP010000002.1": "JACBPP",
    "JACBPP010000003.1": "JACBPP",
    "JACBPP010000005.1": "JACBPP",
    "JACBPP010000011.1": "JACBPP",
    "JACBPP010000012.1": "JACBPP",
    "JACBPP010000013.1": "JACBPP",
    "JACBPP010000014.1": "JACBPP",
    "JACBPP010000015.1": "JACBPP",
    "JACBPP010000016.1": "JACBPP",
    
    # JAJMCQ
    "JAJMCQ010000001.1": "JAJMCQ",
    "JAJMCQ010000002.1": "JAJMCQ",
    "JAJMCQ010000003.1": "JAJMCQ",
    "JAJMCQ010000004.1": "JAJMCQ",
    "JAJMCQ010000005.1": "JAJMCQ",
    "JAJMCQ010000006.1": "JAJMCQ",
    "JAJMCQ010000007.1": "JAJMCQ",
    "JAJMCQ010000008.1": "JAJMCQ",
    "JAJMCQ010000009.1": "JAJMCQ",
    "JAJMCQ010000010.1": "JAJMCQ",
    "JAJMCQ010000011.1": "JAJMCQ",
    "JAJMCQ010000012.1": "JAJMCQ",
    "JAJMCQ010000013.1": "JAJMCQ",
    "JAJMCQ010000014.1": "JAJMCQ",
    "JAJMCQ010000015.1": "JAJMCQ",
    "JAJMCQ010000016.1": "JAJMCQ",
    "JAJMCQ010000017.1": "JAJMCQ",
    "JAJMCQ010000018.1": "JAJMCQ",
    "JAJMCQ010000019.1": "JAJMCQ",
    "JAJMCQ010000020.1": "JAJMCQ",
    
    # MTJM
    "MTJM01000001.1": "MTJM",
    "MTJM01000002.1": "MTJM",
    "MTJM01000003.1": "MTJM",
    "MTJM01000004.1": "MTJM",
    "MTJM01000005.1": "MTJM",
    "MTJM01000006.1": "MTJM",
    "MTJM01000007.1": "MTJM",
    "MTJM01000008.1": "MTJM",
    "MTJM01000009.1": "MTJM",
    "MTJM01000010.1": "MTJM",
    
    # JAKTUL
    "JAKTUL010001563.1": "JAKTUL",
    "JAKTUL010001572.1": "JAKTUL",
    "JAKTUL010001573.1": "JAKTUL",
    "JAKTUL010001574.1": "JAKTUL",
    "JAKTUL010001575.1": "JAKTUL",
    "JAKTUL010001576.1": "JAKTUL",
    "JAKTUL010001577.1": "JAKTUL",
    "JAKTUL010001578.1": "JAKTUL",
    "JAKTUL010001579.1": "JAKTUL",
    "JAKTUL010001580.1": "JAKTUL",
    "JAKTUL010001581.1": "JAKTUL",
    "JAKTUL010001582.1": "JAKTUL",
    "JAKTUL010001583.1": "JAKTUL",
    "JAKTUL010001584.1": "JAKTUL",
    "JAKTUL010001585.1": "JAKTUL",
    "JAKTUL010001586.1": "JAKTUL",
    "JAKTUL010001587.1": "JAKTUL",
    "JAKTUL010001589.1": "JAKTUL",
    "JAKTUL010001590.1": "JAKTUL",
}

# Step 2: Load your annotation file (Prokka/GFF converted to TSV)
annotation_file = r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\GEN.tsv"
df = pd.read_csv(annotation_file, sep="\t")

# Step 3: Map contig IDs to strain
df['strain'] = df['locus_tag'].map(contig_to_strain)

# Step 4: Save new annotated file
output_file = r"C:\Users\Giorgia\dna\dna2vec-legacy\genes_with_strain.tsv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Strain IDs added. Saved to {output_file}")
import pandas as pd

# Load annotation file
df = pd.read_csv(r"C:\Users\Giorgia\Documents\M_pulcherrima_phenotype_prediction\data\GEN.tsv", sep="\t")

# Example contig → strain mapping
contig_to_strain = {
    "CP034462.1": "CP",
    "CP034461.1": "CP",
    "CP034460.1": "CP",
    "CP034459.1": "CP",
    "CP034458.1": "CP",
    "CP034457.1": "CP",
    "CP034456.1": "CP",
    
    # JAQFWG
    "JAQFWG010000001.1": "JAQFWG",
    "JAQFWG010000002.1": "JAQFWG",
    "JAQFWG010000003.1": "JAQFWG",
    "JAQFWG010000006.1": "JAQFWG",
    "JAQFWG010000007.1": "JAQFWG",
    "JAQFWG010000008.1": "JAQFWG",
    "JAQFWG010000009.1": "JAQFWG",
    "JAQFWG010000011.1": "JAQFWG",
    "JAQFWG010000012.1": "JAQFWG",
    "JAQFWG010000013.1": "JAQFWG",
    "JAQFWG010000014.1": "JAQFWG",
    "JAQFWG010000015.1": "JAQFWG",
    "JAQFWG010000016.1": "JAQFWG",
    "JAQFWG010000017.1": "JAQFWG",
    "JAQFWG010000018.1": "JAQFWG",
    "JAQFWG010000019.1": "JAQFWG",
    "JAQFWG010000020.1": "JAQFWG",
    "JAQFWG010000021.1": "JAQFWG",
    "JAQFWG010000022.1": "JAQFWG",
    "JAQFWG010000023.1": "JAQFWG",
    "JAQFWG010000024.1": "JAQFWG",
    
    # JAJMIJ
    "JAJMIJ010000343.1": "JAJMIJ",
    "JAJMIJ010000344.1": "JAJMIJ",
    "JAJMIJ010000345.1": "JAJMIJ",
    "JAJMIJ010000346.1": "JAJMIJ",
    "JAJMIJ010000347.1": "JAJMIJ",
    "JAJMIJ010000348.1": "JAJMIJ",
    "JAJMIJ010000349.1": "JAJMIJ",
    "JAJMIJ010000350.1": "JAJMIJ",
    "JAJMIJ010000351.1": "JAJMIJ",
    "JAJMIJ010000352.1": "JAJMIJ",
    "JAJMIJ010000353.1": "JAJMIJ",
    "JAJMIJ010000354.1": "JAJMIJ",
    "JAJMIJ010000355.1": "JAJMIJ",
    "JAJMIJ010000356.1": "JAJMIJ",
    "JAJMIJ010000357.1": "JAJMIJ",
    "JAJMIJ010000358.1": "JAJMIJ",
    "JAJMIJ010000360.1": "JAJMIJ",
    "JAJMIJ010000361.1": "JAJMIJ",
    "JAJMIJ010000362.1": "JAJMIJ",
    
    # NISE
    "NISE01000342.1": "NISE",
    "NISE01000343.1": "NISE",
    "NISE01000344.1": "NISE",
    "NISE01000345.1": "NISE",
    "NISE01000346.1": "NISE",
    "NISE01000347.1": "NISE",
    "NISE01000348.1": "NISE",
    "NISE01000349.1": "NISE",
    "NISE01000350.1": "NISE",
    "NISE01000351.1": "NISE",
    "NISE01000352.1": "NISE",
    "NISE01000353.1": "NISE",
    "NISE01000354.1": "NISE",
    "NISE01000355.1": "NISE",
    "NISE01000356.1": "NISE",
    "NISE01000357.1": "NISE",
    "NISE01000358.1": "NISE",
    "NISE01000359.1": "NISE",
    "NISE01000360.1": "NISE",
    "NISE01000361.1": "NISE",
    
    # JAKTYT
    "JAKTYT010000001.1": "JAKTYT",
    "JAKTYT010000002.1": "JAKTYT",
    "JAKTYT010000003.1": "JAKTYT",
    "JAKTYT010000005.1": "JAKTYT",
    "JAKTYT010000007.1": "JAKTYT",
    "JAKTYT010000008.1": "JAKTYT",
    "JAKTYT010000009.1": "JAKTYT",
    "JAKTYT010000010.1": "JAKTYT",
    "JAKTYT010000011.1": "JAKTYT",
    "JAKTYT010000012.1": "JAKTYT",
    "JAKTYT010000013.1": "JAKTYT",
    "JAKTYT010000015.1": "JAKTYT",
    "JAKTYT010000017.1": "JAKTYT",
    "JAKTYT010000018.1": "JAKTYT",
    "JAKTYT010000019.1": "JAKTYT",
    "JAKTYT010000020.1": "JAKTYT",
    
    # QBLL
    "QBLL01000001.1": "QBLL",
    "QBLL01000002.1": "QBLL",
    "QBLL01000003.1": "QBLL",
    "QBLL01000004.1": "QBLL",
    "QBLL01000005.1": "QBLL",
    "QBLL01000006.1": "QBLL",
    "QBLL01000007.1": "QBLL",
    "QBLL01000008.1": "QBLL",
    "QBLL01000009.1": "QBLL",
    "QBLL01000010.1": "QBLL",
    "QBLL01000011.1": "QBLL",
    "QBLL01000012.1": "QBLL",
    "QBLL01000013.1": "QBLL",
    "QBLL01000014.1": "QBLL",
    "QBLL01000015.1": "QBLL",
    "QBLL01000016.1": "QBLL",
    "QBLL01000017.1": "QBLL",
    "QBLL01000018.1": "QBLL",
    "QBLL01000019.1": "QBLL",
    "QBLL01000020.1": "QBLL",
    
    # JAJMIQ
    "JAJMIQ010000001.1": "JAJMIQ",
    "JAJMIQ010000002.1": "JAJMIQ",
    "JAJMIQ010000003.1": "JAJMIQ",
    "JAJMIQ010000004.1": "JAJMIQ",
    "JAJMIQ010000005.1": "JAJMIQ",
    "JAJMIQ010000006.1": "JAJMIQ",
    "JAJMIQ010000007.1": "JAJMIQ",
    "JAJMIQ010000008.1": "JAJMIQ",
    "JAJMIQ010000009.1": "JAJMIQ",
    "JAJMIQ010000010.1": "JAJMIQ",
    "JAJMIQ010000011.1": "JAJMIQ",
    "JAJMIQ010000012.1": "JAJMIQ",
    "JAJMIQ010000013.1": "JAJMIQ",
    "JAJMIQ010000014.1": "JAJMIQ",
    "JAJMIQ010000015.1": "JAJMIQ",
    "JAJMIQ010000016.1": "JAJMIQ",
    "JAJMIQ010000017.1": "JAJMIQ",
    "JAJMIQ010000018.1": "JAJMIQ",
    "JAJMIQ010000019.1": "JAJMIQ",
    "JAJMIQ010000020.1": "JAJMIQ",
    
    # JAJMHS
    "JAJMHS010000001.1": "JAJMHS",
    "JAJMHS010000002.1": "JAJMHS",
    "JAJMHS010000003.1": "JAJMHS",
    "JAJMHS010000004.1": "JAJMHS",
    "JAJMHS010000005.1": "JAJMHS",
    "JAJMHS010000006.1": "JAJMHS",
    "JAJMHS010000007.1": "JAJMHS",
    "JAJMHS010000008.1": "JAJMHS",
    "JAJMHS010000009.1": "JAJMHS",
    "JAJMHS010000010.1": "JAJMHS",
    "JAJMHS010000011.1": "JAJMHS",
    "JAJMHS010000012.1": "JAJMHS",
    "JAJMHS010000013.1": "JAJMHS",
    "JAJMHS010000014.1": "JAJMHS",
    "JAJMHS010000015.1": "JAJMHS",
    "JAJMHS010000016.1": "JAJMHS",
    "JAJMHS010000017.1": "JAJMHS",
    "JAJMHS010000018.1": "JAJMHS",
    "JAJMHS010000019.1": "JAJMHS",
    "JAJMHS010000020.1": "JAJMHS",
    
    # JAKTYM
    "JAKTYM010000001.1": "JAKTYM",
    "JAKTYM010000002.1": "JAKTYM",
    "JAKTYM010000003.1": "JAKTYM",
    "JAKTYM010000004.1": "JAKTYM",
    "JAKTYM010000005.1": "JAKTYM",
    "JAKTYM010000006.1": "JAKTYM",
    "JAKTYM010000007.1": "JAKTYM",
    "JAKTYM010000008.1": "JAKTYM",
    "JAKTYM010000009.1": "JAKTYM",
    "JAKTYM010000010.1": "JAKTYM",
    "JAKTYM010000011.1": "JAKTYM",
    "JAKTYM010000012.1": "JAKTYM",
    "JAKTYM010000013.1": "JAKTYM",
    "JAKTYM010000014.1": "JAKTYM",
    "JAKTYM010000015.1": "JAKTYM",
    "JAKTYM010000016.1": "JAKTYM",
    "JAKTYM010000017.1": "JAKTYM",
    "JAKTYM010000018.1": "JAKTYM",
    "JAKTYM010000019.1": "JAKTYM",
    "JAKTYM010000020.1": "JAKTYM",
    
    # JACBPP
    "JACBPP010000001.1": "JACBPP",
    "JACBPP010000002.1": "JACBPP",
    "JACBPP010000003.1": "JACBPP",
    "JACBPP010000005.1": "JACBPP",
    "JACBPP010000011.1": "JACBPP",
    "JACBPP010000012.1": "JACBPP",
    "JACBPP010000013.1": "JACBPP",
    "JACBPP010000014.1": "JACBPP",
    "JACBPP010000015.1": "JACBPP",
    "JACBPP010000016.1": "JACBPP",
    
    # JAJMCQ
    "JAJMCQ010000001.1": "JAJMCQ",
    "JAJMCQ010000002.1": "JAJMCQ",
    "JAJMCQ010000003.1": "JAJMCQ",
    "JAJMCQ010000004.1": "JAJMCQ",
    "JAJMCQ010000005.1": "JAJMCQ",
    "JAJMCQ010000006.1": "JAJMCQ",
    "JAJMCQ010000007.1": "JAJMCQ",
    "JAJMCQ010000008.1": "JAJMCQ",
    "JAJMCQ010000009.1": "JAJMCQ",
    "JAJMCQ010000010.1": "JAJMCQ",
    "JAJMCQ010000011.1": "JAJMCQ",
    "JAJMCQ010000012.1": "JAJMCQ",
    "JAJMCQ010000013.1": "JAJMCQ",
    "JAJMCQ010000014.1": "JAJMCQ",
    "JAJMCQ010000015.1": "JAJMCQ",
    "JAJMCQ010000016.1": "JAJMCQ",
    "JAJMCQ010000017.1": "JAJMCQ",
    "JAJMCQ010000018.1": "JAJMCQ",
    "JAJMCQ010000019.1": "JAJMCQ",
    "JAJMCQ010000020.1": "JAJMCQ",
    
    # MTJM
    "MTJM01000001.1": "MTJM",
    "MTJM01000002.1": "MTJM",
    "MTJM01000003.1": "MTJM",
    "MTJM01000004.1": "MTJM",
    "MTJM01000005.1": "MTJM",
    "MTJM01000006.1": "MTJM",
    "MTJM01000007.1": "MTJM",
    "MTJM01000008.1": "MTJM",
    "MTJM01000009.1": "MTJM",
    "MTJM01000010.1": "MTJM",
    
    # JAKTUL
    "JAKTUL010001563.1": "JAKTUL",
    "JAKTUL010001572.1": "JAKTUL",
    "JAKTUL010001573.1": "JAKTUL",
    "JAKTUL010001574.1": "JAKTUL",
    "JAKTUL010001575.1": "JAKTUL",
    "JAKTUL010001576.1": "JAKTUL",
    "JAKTUL010001577.1": "JAKTUL",
    "JAKTUL010001578.1": "JAKTUL",
    "JAKTUL010001579.1": "JAKTUL",
    "JAKTUL010001580.1": "JAKTUL",
    "JAKTUL010001581.1": "JAKTUL",
    "JAKTUL010001582.1": "JAKTUL",
    "JAKTUL010001583.1": "JAKTUL",
    "JAKTUL010001584.1": "JAKTUL",
    "JAKTUL010001585.1": "JAKTUL",
    "JAKTUL010001586.1": "JAKTUL",
    "JAKTUL010001587.1": "JAKTUL",
    "JAKTUL010001589.1": "JAKTUL",
    "JAKTUL010001590.1": "JAKTUL",
}

# Add contig column if your locus_tag contains contig info
# Example: locus_tag = ICLPPAIG_00001 might belong to CP034462.1
# If locus_tag mapping exists:
df["contig"] = df["locus_tag"].apply(lambda x: "CP034462.1")  # replace with actual logic

# Map contig → strain
df["strain"] = df["contig"].map(contig_to_strain)

# Check
print(df.head())
print(df["strain"].unique())

