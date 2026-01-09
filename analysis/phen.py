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
import csv

# List of your strains and their corresponding GenBank files
strain_files = {
    "VFXK": "VFXK.gbk",
    "CP": "CP.gbk",
    "JACBPP": "JACBPP.gbk",
    "ANFW": "ANFW.gbk",
    "JAQFWG": "JAQFWG.gbk",
    "JAJMHS": "JAJMHS.gbk",
    "MTJM": "MTJM.gbk",
    "JAKTYM": "JAKTYM.gbk",
    "JAJMCQ": "JAJMCQ.gbk",
    "JAJMIJ": "JAJMIJ.gbk",
    "JAJMIQ": "JAJMIQ.gbk",
    "QBLL": "QBLL.gbk",
    "JAKTYT": "JAKTYT.gbk",
    "NISE": "NISE.gbk",
    "JAKTUL": "JAKTUL.gbk"
}

# Load locus_tag list (from your CSV or TSV)
locus_tags_file = "genes.csv"  # file with column: locus_tag
with open(locus_tags_file, newline='') as f:
    reader = csv.DictReader(f)
    locus_tags_list = [row['locus_tag'] for row in reader]

# Prepare output
output_file = "genes_mapped.csv"
with open(output_file, 'w', newline='') as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow(["Strain", "Contig", "Locus_Tag", "Product", "Functional_Category"])
    
    # Iterate over each strain
    for strain, gbk_file in strain_files.items():
        for record in SeqIO.parse(gbk_file, "genbank"):
            contig_id = record.id
            for feature in record.features:
                if feature.type == "gene" or feature.type == "CDS":
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    product = feature.qualifiers.get("product", [""])[0]
                    functional_category = feature.qualifiers.get("function", ["Other"])[0]
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
df = pd.read_csv(r"C:\Users\Giorgia\dna\dna2vec-legacy\representative_genomes\GEN_annotation\GEN.tsv", sep="\t")

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
                    
                    if locus_tag in locus_tags_list:
                        writer.writerow([strain, contig_id, locus_tag, product, functional_category])

print(f"Mapping complete! Output saved to {output_file}")
