import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# =================================================
# 1. Load CSV
# =================================================
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

# =================================================
# 2. Low variance noise
# =================================================
low_variance_traits = [
    "H2S","Esterase","SO2_Resistance","Sugar_Utilization",
    "Lipase","Growth_Rate","Glucosidase","Gluconic_Acid_Utilization"
]

np.random.seed(42)
for t in low_variance_traits:
    if df[t].std() < 1:
        df[t] += np.random.normal(0, 0.2, size=len(df))

# =================================================
# 3. Assign strain types
# =================================================
def assign_type(strain):
    if "CBS" in strain or "NRRL" in strain:
        return "Type strain"
    elif "AP" in strain or "KIOM" in strain:
        return "Environmental"
    else:
        return "Clinical"

df["Strain_Type"] = df["Strain"].apply(assign_type)

# =================================================
# 4. Standardize traits (base year = 2025)
# =================================================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[traits])

df_scaled = pd.DataFrame(X_scaled, columns=traits)
df_scaled["Strain"] = df["Strain"]
df_scaled["Strain_Type"] = df["Strain_Type"]

# =================================================
# 5. Simulate yearly drift (2015–2035)
# =================================================
years = np.arange(2015, 2036)
base_year = 2025
change_fraction = 0.05

yearly_data = {}
yearly_data[base_year] = df_scaled.copy()

# Past years
for year in range(base_year-1, years.min()-1, -1):
    prev_df = yearly_data[year+1].copy()
    for t in traits:
        std_t = df_scaled[t].std()
        prev_df[t] -= np.random.normal(0, std_t*change_fraction, size=len(prev_df))
    yearly_data[year] = prev_df

# Future years
for year in range(base_year+1, years.max()+1):
    prev_df = yearly_data[year-1].copy()
    for t in traits:
        std_t = df_scaled[t].std()
        prev_df[t] += np.random.normal(0, std_t*change_fraction, size=len(prev_df))
    yearly_data[year] = prev_df

# =================================================
# 6. PCA on all years
# =================================================
pca = PCA(n_components=2, random_state=42)
all_years_matrix = pd.concat([yearly_data[y][traits] for y in years], ignore_index=True)
pca.fit(all_years_matrix)

for y in years:
    yearly_data[y][["PC1","PC2"]] = pca.transform(yearly_data[y][traits])

# =================================================
# 7. Animation setup
# =================================================
fig, ax = plt.subplots(figsize=(12,8))
palette = {
    "Type strain":"tab:blue",
    "Environmental":"tab:green",
    "Clinical":"tab:red"
}

points = {strain: ax.plot([], [], 'o', color=palette[stype], markersize=8)[0]
          for strain, stype in zip(df_scaled["Strain"], df_scaled["Strain_Type"])}

labels = {strain: ax.text(0,0,strain,fontsize=7) for strain in df_scaled["Strain"]}

frames_per_year = 10
total_frames = (len(years)-1)*frames_per_year

# Year stamp (timer) text
year_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                    fontsize=16, fontweight='bold', color='black')

# =================================================
# 8. Init function
# =================================================
def init():
    pcs = pd.concat([yearly_data[y][["PC1","PC2"]] for y in years])
    ax.set_xlim(pcs["PC1"].min()-1, pcs["PC1"].max()+1)
    ax.set_ylim(pcs["PC2"].min()-1, pcs["PC2"].max()+1)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Yearly PCA Trajectories of Yeast Strains (2015–2035)")
    return list(points.values()) + list(labels.values()) + [year_text]

# =================================================
# 9. Update function
# =================================================
def update(frame):
    year_idx = frame // frames_per_year
    intra = (frame % frames_per_year) / frames_per_year

    y0 = years[year_idx]
    y1 = years[min(year_idx+1,len(years)-1)]
    current_year = y0 + intra

    for strain in df_scaled["Strain"]:
        p0 = yearly_data[y0].loc[yearly_data[y0]["Strain"]==strain, ["PC1","PC2"]].values[0]
        p1 = yearly_data[y1].loc[yearly_data[y1]["Strain"]==strain, ["PC1","PC2"]].values[0]

        x = p0[0] + intra*(p1[0]-p0[0])
        y = p0[1] + intra*(p1[1]-p0[1])

        points[strain].set_data(x, y)
        labels[strain].set_position((x+0.04, y+0.04))

    # Update year stamp
    year_text.set_text(f"Year: {current_year:.1f}")

    return list(points.values()) + list(labels.values()) + [year_text]

# =================================================
# 10. Run animation
# =================================================
ani = FuncAnimation(fig, update, frames=total_frames, init_func=init,
                    blit=True, repeat=True)

plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text
from matplotlib import patheffects

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

numeric_traits = [col for col in traits if df[col].dtype in ["float64", "int64"]]
categorical_traits = [col for col in traits if col not in numeric_traits]

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

    # PCA
    if len(numeric_in_group) > 0:
        X_scaled = StandardScaler().fit_transform(df[numeric_in_group])
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(X_scaled)
        df["PC1"] = pcs[:, 0]
        df["PC2"] = pcs[:, 1]
    else:
        df["PC1"], df["PC2"] = 0, 0

    # Add small jitter to avoid straight line collapse
    jitter_strength = 0.03
    df["PC1_plot"] = df["PC1"] + np.random.uniform(-jitter_strength, jitter_strength, size=df.shape[0])
    df["PC2_plot"] = df["PC2"] + np.random.uniform(-jitter_strength, jitter_strength, size=df.shape[0])

    n = len(group_traits)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for idx, trait in enumerate(group_traits):
        ax = axes[idx]
        marker_size = 60

        if trait in numeric_traits:
            scatter = ax.scatter(
                df["PC1_plot"], df["PC2_plot"], c=df[trait], cmap="viridis",
                s=marker_size, edgecolors='k', linewidth=0.3
            )
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.ax.tick_params(labelsize=6)
            cbar.set_label(trait, fontsize=8)
        else:
            sns.scatterplot(
                data=df,
                x="PC1_plot",
                y="PC2_plot",
                hue=trait,
                palette="tab10",
                s=marker_size,
                ax=ax,
                legend=False,
                edgecolor='k',
                linewidth=0.3
            )

        # ----------------------------
        # Strain labels with strong halo
        # ----------------------------
        texts = []
        for i in range(len(df)):
            texts.append(
                ax.text(
                    df["PC1_plot"][i], df["PC2_plot"][i], df["Strain"][i],
                    fontsize=6,  # larger font for visibility
                    fontweight='bold',
                    alpha=0.85,
                    path_effects=[
                        patheffects.Stroke(linewidth=2, foreground="white"),
                        patheffects.Normal()
                    ]
                )
            )

        adjust_text(
            texts,
            ax=ax,
            only_move={'points':'y','texts':'y'},
            force_text=0.7,
            force_points=0.5,
            expand_text=(1.1, 2.0),
            arrowprops=dict(arrowstyle="-", color='gray', lw=0.2)
        )

        ax.set_title(trait, fontsize=12)
        ax.set_xlabel("PC1", fontsize=8)
        ax.set_ylabel("PC2", fontsize=8)
        ax.tick_params(axis="both", labelsize=6)
        ax.grid(True, linestyle='--', alpha=0.3)

    plt.suptitle(f"PCA Group: {group_name}", fontsize=15)
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import networkx as nx

# -----------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

# -----------------------------
# Traits and Origin
# -----------------------------
traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

def assign_origin(strain):
    if "CBS" in strain or "NRRL" in strain:
        return "Type"
    elif "AP" in strain or "KIOM" in strain:
        return "Environmental"
    else:
        return "Clinical"

df["Origin"] = df["Strain"].apply(assign_origin)

# -----------------------------
# Standardize
# -----------------------------
X_scaled = StandardScaler().fit_transform(df[traits])
df_scaled = pd.DataFrame(X_scaled, columns=traits)
df_scaled["Strain"] = df["Strain"]
df_scaled["Origin"] = df["Origin"]

# -----------------------------
# PCA
# -----------------------------
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(X_scaled)
df_scaled["PC1"], df_scaled["PC2"] = pca_coords.T

# -----------------------------
# Correlation Network
# -----------------------------
corr_matrix = df_scaled[traits].corr()
edges = [(i,j) for i in traits for j in traits if i!=j and abs(corr_matrix.loc[i,j])>=0.5]

# Define node groups
node_colors = []
for t in traits:
    if t in ["Sugar_Utilization","Growth_Rate","Gluconic_Acid_Utilization","Glucosidase"]:
        node_colors.append("green")  # Metabolic core
    elif t in ["Ethanol_Tolerance","Copper_Resistance","SO2_Resistance"]:
        node_colors.append("orange") # Stress traits
    else:
        node_colors.append("blue")   # Other (H2S, etc.)

G = nx.Graph()
G.add_nodes_from(traits)
G.add_edges_from(edges)

# -----------------------------
# Plotting
# -----------------------------
fig, axes = plt.subplots(1, 2, figsize=(16,7))

# PCA plot
sns.scatterplot(
    data=df_scaled, x="PC1", y="PC2", hue="Origin",
    palette={"Type":"red","Environmental":"green","Clinical":"blue"},
    ax=axes[0], s=100, edgecolor="k"
)
axes[0].set_title("PCA of M. pulcherrima Traits", fontsize=14)
axes[0].legend(title="Origin")
axes[0].grid(True)

# Correlation network
pos = nx.circular_layout(G)
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=800, ax=axes[1])
# Edge style based on correlation magnitude
for i,j in G.edges():
    weight = abs(corr_matrix.loc[i,j])
    style = 'solid' if weight>=0.7 else 'dashed'
    nx.draw_networkx_edges(G, pos, edgelist=[(i,j)], style=style, width=2, ax=axes[1])
nx.draw_networkx_labels(G, pos, font_size=10, ax=axes[1])
axes[1].set_title("Trait Correlation Network", fontsize=14)
axes[1].axis('off')

plt.tight_layout()
plt.show()

