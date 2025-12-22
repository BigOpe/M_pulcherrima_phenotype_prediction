import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text
from matplotlib import patheffects   # <-- FIXED

# ------------------------------------
# Load data
# ------------------------------------
df = pd.read_csv(r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv")

traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

# ------------------------------------
# Numeric vs categorical traits
# ------------------------------------
numeric_traits = [col for col in traits if df[col].dtype in ["float64", "int64"]]
categorical_traits = [col for col in traits if col not in numeric_traits]

# ------------------------------------
# Functional groups
# ------------------------------------
trait_groups = {
    "Stress Tolerance": ["Ethanol_Tolerance", "Copper_Resistance", "SO2_Resistance"],
    "Metabolism / Sugar Utilization": [
        "Sugar_Utilization", "Glucosidase", "Gluconic_Acid_Utilization", "Esterase"
    ],
    "Growth / Enzymatic Activity": [
        "Growth_Rate", "Protease", "Lipase", 
        "Killer_Activity", "Pigment_Intensity", "H2S"
    ]
}


# ------------------------------------
# PCA & Plotting
# ------------------------------------
for group_name, group_traits in trait_groups.items():

    numeric_in_group = [t for t in group_traits if t in numeric_traits]

    # PCA on numeric traits only
    if len(numeric_in_group) > 0:
        X_scaled = StandardScaler().fit_transform(df[numeric_in_group])
        pca = PCA(n_components=2)
        df["PC1"], df["PC2"] = pca.fit_transform(X_scaled).T
    else:
        df["PC1"], df["PC2"] = 0, 0  # unlikely, but safe fallback

    # Create subplot grid based on number of traits
    n = len(group_traits)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))

    if n == 1:
        axes = [axes]

    # ------------------------------------
    # Plot each trait
    # ------------------------------------
    for idx, trait in enumerate(group_traits):
        ax = axes[idx]

        marker_size = 60

        # Numeric trait color map
        if trait in numeric_traits:
            scatter = ax.scatter(
                df["PC1"], df["PC2"], c=df[trait], cmap="viridis", s=marker_size
            )
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.ax.tick_params(labelsize=6)
            cbar.set_label(trait, fontsize=8)

        # Categorical trait color map
        else:
            sns.scatterplot(
                data=df,
                x="PC1",
                y="PC2",
                hue=trait,
                palette="tab10",
                s=marker_size,
                ax=ax,
                legend=False
            )

        # ----------------------------
        # Strain labels with halo outline
        # ----------------------------
        texts = []
        for i in range(len(df)):
            texts.append(
                ax.text(
                    df["PC1"][i], df["PC2"][i], df["Strain"][i],
                    fontsize=4, alpha=0.65,
                    path_effects=[
                        patheffects.Stroke(linewidth=1, foreground="white"),
                        patheffects.Normal()
                    ]
                )
            )

        adjust_text(
            texts,
            ax=ax,
            only_move={'points': 'y', 'texts': 'y'},
            force_points=0.1,
            force_text=0.1
        )

        ax.set_title(trait, fontsize=11)
        ax.set_xlabel("PC1", fontsize=8)
        ax.set_ylabel("PC2", fontsize=8)
        ax.tick_params(axis="both", labelsize=7)

    plt.suptitle(f"PCA Group: {group_name}", fontsize=15)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.show()

# ------------------------------------
# DENDROGRAM (Hierarchical Clustering)
# ------------------------------------
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

# Use all numeric traits for clustering
numeric_data = df[numeric_traits].copy()

# Standardize
scaler = StandardScaler()
scaled = scaler.fit_transform(numeric_data)

# Compute hierarchical clustering
linkage_matrix = linkage(scaled, method='ward')

# Plot dendrogram
plt.figure(figsize=(12, 6))
dendrogram(
    linkage_matrix,
    labels=df["Strain"].values,
    leaf_rotation=90,
    leaf_font_size=8,
)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Strain")
plt.ylabel("Distance")
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text
from matplotlib import patheffects
from matplotlib.patches import Patch
from scipy.spatial import ConvexHull

# ------------------------------------
# Load data
# ------------------------------------
df = pd.read_csv(r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv")

traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

# ------------------------------------
# Assign strain types
# ------------------------------------
def assign_type(strain_name):
    if "CBS" in strain_name or "NRRL" in strain_name:
        return "Type strain"
    elif "AP" in strain_name or "KIOM" in strain_name:
        return "Environmental"
    else:
        return "Clinical"

df["Strain_Type"] = df["Strain"].apply(assign_type)

# ------------------------------------
# Trait groups
# ------------------------------------
trait_groups = {
    "Stress Tolerance": ["Ethanol_Tolerance", "Copper_Resistance", "SO2_Resistance"],
    "Metabolism / Sugar Utilization": [
        "Sugar_Utilization", "Glucosidase", "Gluconic_Acid_Utilization", "Esterase"
    ],
    "Growth / Enzymatic Activity": [
        "Growth_Rate", "Protease", "Lipase", 
        "Killer_Activity", "Pigment_Intensity", "H2S"
    ]
}

# ------------------------------------
# Build strain-to-color mapping
# ------------------------------------
clinical_strains = df[df["Strain_Type"]=="Clinical"]["Strain"].unique()
env_strains      = df[df["Strain_Type"]=="Environmental"]["Strain"].unique()
type_strains     = df[df["Strain_Type"]=="Type strain"]["Strain"].unique()

clinical_palette = sns.color_palette("Reds", n_colors=len(clinical_strains))
env_palette      = sns.color_palette("Greens", n_colors=len(env_strains))
type_palette     = sns.color_palette("Blues", n_colors=len(type_strains))

strain_colors = {}
for strain, color in zip(clinical_strains, clinical_palette):
    strain_colors[strain] = color
for strain, color in zip(env_strains, env_palette):
    strain_colors[strain] = color
for strain, color in zip(type_strains, type_palette):
    strain_colors[strain] = color

# ------------------------------------
# Helper: draw convex hull for each group
# ------------------------------------
def draw_hull(ax, data, group, color):
    points = data[df["Strain_Type"]==group][["PC1","PC2"]].values
    if len(points) >= 3:  # need at least 3 points for a hull
        hull = ConvexHull(points)
        hull_points = points[hull.vertices]
        ax.fill(hull_points[:,0], hull_points[:,1], alpha=0.15, color=color, label=f"{group} hull")

# ------------------------------------
# PCA plots for each trait group
# ------------------------------------
for group_name, group_traits in trait_groups.items():
    # numeric traits only
    numeric_in_group = [t for t in group_traits if df[t].dtype in ["float64","int64"]]
    if len(numeric_in_group) > 0:
        X_scaled = StandardScaler().fit_transform(df[numeric_in_group])
        pca = PCA(n_components=2)
        df["PC1"], df["PC2"] = pca.fit_transform(X_scaled).T
    else:
        df["PC1"], df["PC2"] = 0, 0

    plt.figure(figsize=(8,6))
    ax = plt.gca()

    # scatter each strain with its unique color
    for i in range(len(df)):
        ax.scatter(
            df["PC1"][i], df["PC2"][i],
            color=strain_colors[df["Strain"][i]],
            s=70, edgecolor="black", linewidth=0.5
        )
        ax.text(
            df["PC1"][i], df["PC2"][i], df["Strain"][i],
            fontsize=6, alpha=0.7,
            path_effects=[
                patheffects.Stroke(linewidth=1, foreground="white"),
                patheffects.Normal()
            ]
        )

    # convex hulls for each type
    draw_hull(ax, df, "Clinical", "red")
    draw_hull(ax, df, "Environmental", "green")
    draw_hull(ax, df, "Type strain", "blue")

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"PCA Group: {group_name}")

    # legend for strain types
    handles = [
        Patch(facecolor="red", label="Clinical"),
        Patch(facecolor="green", label="Environmental"),
        Patch(facecolor="blue", label="Type strain")
    ]
    ax.legend(handles=handles, title="Strain Type", bbox_to_anchor=(1.05,1), loc="upper left")

    plt.tight_layout()
    plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text
from matplotlib import patheffects
from matplotlib.patches import Patch

# ------------------------------------
# Load data
# ------------------------------------
df = pd.read_csv(r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv")

traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

# ------------------------------------
# Assign strain types
# ------------------------------------
def assign_type(strain_name):
    if "CBS" in strain_name or "NRRL" in strain_name:
        return "Type strain"
    elif "AP" in strain_name or "KIOM" in strain_name:
        return "Environmental"
    else:
        return "Clinical"

df["Strain_Type"] = df["Strain"].apply(assign_type)

# ------------------------------------
# Functional groups
# ------------------------------------
trait_groups = {
    "Stress Tolerance": ["Ethanol_Tolerance", "Copper_Resistance", "SO2_Resistance"],
    "Metabolism / Sugar Utilization": [
        "Sugar_Utilization", "Glucosidase", "Gluconic_Acid_Utilization", "Esterase"
    ],
    "Growth / Enzymatic Activity": [
        "Growth_Rate", "Protease", "Lipase", 
        "Killer_Activity", "Pigment_Intensity", "H2S"
    ]
}

# ------------------------------------
# Color palette for strain types
# ------------------------------------
type_colors = {"Clinical":"red","Environmental":"green","Type strain":"blue"}

# ------------------------------------
# PCA & Plotting
# ------------------------------------
# ------------------------------------
# PCA & Plotting (IMPROVED LABEL VISIBILITY)
# ------------------------------------
for group_name, group_traits in trait_groups.items():

    numeric_in_group = [t for t in group_traits if df[t].dtype in ["float64", "int64"]]

    if len(numeric_in_group) > 0:
        X_scaled = StandardScaler().fit_transform(df[numeric_in_group])
        pca = PCA(n_components=2)
        df["PC1"], df["PC2"] = pca.fit_transform(X_scaled).T
    else:
        df["PC1"], df["PC2"] = 0, 0

    # Add slight jitter to avoid perfect overlap
    jitter = 0.02
    df["PC1_j"] = df["PC1"] + np.random.normal(0, jitter, size=len(df))
    df["PC2_j"] = df["PC2"] + np.random.normal(0, jitter, size=len(df))

    plt.figure(figsize=(8, 7))
    ax = plt.gca()

    texts = []

    for i in range(len(df)):
        ax.scatter(
            df["PC1_j"][i], df["PC2_j"][i],
            color=type_colors[df["Strain_Type"][i]],
            s=90,
            edgecolor="black",
            linewidth=0.6,
            zorder=3
        )

        texts.append(
            ax.text(
                df["PC1_j"][i], df["PC2_j"][i],
                df["Strain"][i],
                fontsize=8.5,          # 🔼 BIGGER FONT
                alpha=0.9,
                zorder=4,
                path_effects=[
                    patheffects.Stroke(linewidth=1.8, foreground="white"),
                    patheffects.Normal()
                ]
            )
        )

    # 🔥 Automatic label repulsion
    adjust_text(
        texts,
        ax=ax,
        expand_points=(1.4, 1.6),
        expand_text=(1.3, 1.5),
        arrowprops=dict(arrowstyle="-", lw=0.4, color="gray"),
        force_points=0.6,
        force_text=0.8
    )

    # Legend
    handles = [Patch(facecolor=color, label=stype) for stype, color in type_colors.items()]
    ax.legend(
        handles=handles,
        title="Origin",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        frameon=False
    )

    ax.set_title(f"PCA Group: {group_name}", fontsize=14, pad=12)
    ax.set_xlabel("PC1", fontsize=11)
    ax.set_ylabel("PC2", fontsize=11)
    ax.tick_params(axis="both", labelsize=9)

    ax.grid(alpha=0.2, linestyle="--")
    plt.tight_layout()
    plt.show()
