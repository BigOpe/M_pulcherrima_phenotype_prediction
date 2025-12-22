import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib.patches import Patch

# -----------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

# -----------------------------
# Traits
# -----------------------------
traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

# -----------------------------
# Assign strain types
# -----------------------------
def assign_type(strain_name):
    if "CBS" in strain_name or "NRRL" in strain_name:
        return "Type strain"
    elif "AP" in strain_name or "KIOM" in strain_name:
        return "Environmental"
    else:
        return "Clinical"

df["Strain_Type"] = df["Strain"].apply(assign_type)

# -----------------------------
# Standardize traits
# -----------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[traits])
df_scaled = pd.DataFrame(X_scaled, columns=traits)
df_scaled["Strain"] = df["Strain"]
df_scaled["Strain_Type"] = df["Strain_Type"]

# -----------------------------
# Hierarchical clustering
# -----------------------------
linkage_matrix = linkage(X_scaled, method='ward')

# -----------------------------
# Dendrogram with colored strain types and legend
# -----------------------------
type_colors = {"Clinical":"red","Environmental":"green","Type strain":"blue"}

plt.figure(figsize=(14,7))
dendro = dendrogram(
    linkage_matrix,
    labels=df_scaled["Strain"].values,
    leaf_rotation=90,
    leaf_font_size=10,
)

ax = plt.gca()
xlbls = ax.get_xmajorticklabels()
for lbl in xlbls:
    strain_name = lbl.get_text()
    strain_type = df_scaled.loc[df_scaled["Strain"]==strain_name, "Strain_Type"].values[0]
    lbl.set_color(type_colors[strain_type])

# Legend for strain types
handles = [Patch(facecolor=color, label=stype) for stype, color in type_colors.items()]
plt.legend(handles=handles, title="Origin", bbox_to_anchor=(1.05,1), loc="upper left")

plt.title("Hierarchical Clustering Dendrogram (Strain Type Colored)")
plt.xlabel("Strain")
plt.ylabel("Distance")
plt.tight_layout()
plt.show()

# -----------------------------
# -----------------------------
# Heatmap to show trait patterns per strain
# -----------------------------
row_colors = df_scaled["Strain_Type"].map(type_colors)

g = sns.clustermap(
    df_scaled[traits],
    row_linkage=linkage_matrix,
    col_cluster=False,
    row_colors=row_colors,
    cmap="vlag",
    figsize=(10,12),
    yticklabels=df_scaled["Strain"]
)

# 🔽 Make strain (Y-axis) labels smaller
g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_yticklabels(),
    fontsize=7
)

# Optional: make trait (X-axis) labels smaller too
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    fontsize=7,
    rotation=90,
    ha="right"
)

plt.suptitle(
    "Heatmap of Strain Traits (Grouped by Hierarchical Clustering)",
    y=1.02
)

# Legend for strain types
for stype, color in type_colors.items():
    plt.scatter([], [], color=color, label=stype)
plt.legend(title="Origin", bbox_to_anchor=(1.05,1), loc="upper left")

plt.show()

