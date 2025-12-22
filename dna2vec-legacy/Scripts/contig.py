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
                    
                    if locus_tag in locus_tags_list:
                        writer.writerow([strain, contig_id, locus_tag, product, functional_category])

print(f"Mapping complete! Output saved to {output_file}")
