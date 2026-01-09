import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from adjustText import adjust_text

# -------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

# -----------------------------
# Traits and low-variance imputation
# -----------------------------
traits = ["Ethanol_Tolerance","Killer_Activity","Pigment_Intensity","Copper_Resistance",
          "Sugar_Utilization","H2S","Protease","Lipase","Growth_Rate","Glucosidase",
          "SO2_Resistance","Esterase","Gluconic_Acid_Utilization"]

# Add small random noise to low-variance traits
low_variance_traits = ["H2S","Esterase","SO2_Resistance","Sugar_Utilization","Lipase",
                       "Growth_Rate","Glucosidase","Gluconic_Acid_Utilization"]
np.random.seed(42)
for t in low_variance_traits:
    if df[t].std() < 1:
        df[t] += np.random.normal(0, 0.2, size=len(df))  # small noise

# -----------------------------
# Standardize traits
# -----------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(df[traits])
df_scaled = pd.DataFrame(X_scaled, columns=traits)
df_scaled["Strain"] = df["Strain"]

# -----------------------------
# Add strain type categorical variable
# -----------------------------
def assign_type(strain_name):
    if "CBS" in strain_name or "NRRL" in strain_name:
        return "Type strain"
    elif "AP" in strain_name or "KIOM" in strain_name:
        return "Environmental"
    else:
        return "Clinical"

df_scaled["Strain_Type"] = df_scaled["Strain"].apply(assign_type)

# -----------------------------
# PCA for visualization
# -----------------------------
pca = PCA(n_components=2, random_state=42)
df_scaled[["PC1","PC2"]] = pca.fit_transform(df_scaled[traits])
print("Explained variance ratio (PCA):", pca.explained_variance_ratio_)

plt.figure(figsize=(10,7))

ax = sns.scatterplot(
    data=df_scaled,
    x="PC1",
    y="PC2",
    hue="Strain_Type",
    palette="tab10",
    s=100
)

# -----------------------------
# Add strain labels
# -----------------------------
texts = []
for _, row in df_scaled.iterrows():
    texts.append(
        plt.text(
            row["PC1"],
            row["PC2"],
            row["Strain"],
            fontsize=8
        )
    )

# Automatically adjust text to avoid overlaps
adjust_text(
    texts,
    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
)

plt.title("PCA of Yeast Strains Colored by Type")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
plt.tight_layout()
plt.show()


# -----------------------------
# Regression analysis and plotting
# -----------------------------
# -----------------------------
# Regression analysis and plotting
# -----------------------------
results = []

for trait in traits:
    y = df_scaled[trait].values
    X = df_scaled.drop(traits + ["Strain"], axis=1)  # keep categorical columns only
    # One-hot encode categorical
    X = pd.get_dummies(X, columns=["Strain_Type"], drop_first=True)
    
    # ---- PLS Regression ----
    pls = PLSRegression(n_components=2)
    y_pred_pls = cross_val_predict(pls, X, y, cv=5)
    r2_pls = r2_score(y, y_pred_pls)
    r_pls, p_pls = pearsonr(y, y_pred_pls)
    
    # ---- Random Forest Regression ----
    rf = RandomForestRegressor(n_estimators=500, random_state=42)
    rf.fit(X, y)
    y_pred_rf = rf.predict(X)
    r2_rf = r2_score(y, y_pred_rf)
    r_rf, p_rf = pearsonr(y, y_pred_rf)
    
    results.append({
        "Trait": trait,
        "R2_PLS": r2_pls,
        "r_PLS": r_pls,
        "p_PLS": p_pls,
        "R2_RF": r2_rf,
        "r_RF": r_rf,
        "p_RF": p_rf
    })
    
    # ---- Regression plot (Observed vs Predicted) ----
    plt.figure(figsize=(8,6))
    sns.scatterplot(x=y, y=y_pred_pls, hue=df_scaled["Strain_Type"], palette="tab10", s=80)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], "r--")
    plt.xlabel("Observed")
    plt.ylabel("Predicted")
    
    # Format p-value in scientific notation with minimum threshold
    min_p = 1e-16
    p_display = max(p_pls, min_p)
    p_text = f"p={p_display:.3e}"
    
    plt.title(f"{trait} Regression (PLS)\nR²={r2_pls:.2f}, r={r_pls:.2f}, {p_text}")
    
    # Add strain labels to points
    texts = [plt.text(xo, yp, s, fontsize=7) for xo, yp, s in zip(y, y_pred_pls, df_scaled["Strain"])]
    adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))
    
    plt.tight_layout()
    plt.show()

# -----------------------------
# Results DataFrame
# -----------------------------
df_results = pd.DataFrame(results)
print(df_results)

# -----------------------------
# Save improved CSV
# -----------------------------
output_csv = r"C:\Users\Giorgia\dna\dna2vec-legacy\phen_scaled.csv"
df_scaled.to_csv(output_csv, index=False)
print(f"Improved CSV saved at {output_csv}")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# -----------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\phen_scaled.csv"
df = pd.read_csv(csv_path)

traits = ["Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
          "Copper_Resistance","Sugar_Utilization","H2S","Protease",
          "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
          "Esterase","Gluconic_Acid_Utilization"]

# -----------------------------
# Simulate 10-year future projection
# -----------------------------
years = 10
change_fraction = 0.05  # 5% of trait SD per year

future_df = df.copy()
np.random.seed(42)  # reproducibility

for t in traits:
    std_t = df[t].std()
    yearly_change = std_t * change_fraction
    future_df[t] = df[t] + np.random.normal(0, yearly_change, size=len(df)) * years

# -----------------------------
# Standardize traits
# -----------------------------
scaler = StandardScaler()
# Combine current and future for a shared PCA space
combined_df = pd.concat([df, future_df], ignore_index=True)
X_scaled = scaler.fit_transform(combined_df[traits])
combined_df_scaled = pd.DataFrame(X_scaled, columns=traits)
combined_df_scaled["Strain"] = pd.concat([df["Strain"], future_df["Strain"]], ignore_index=True)
combined_df_scaled["Strain_Type"] = pd.concat([df["Strain_Type"], future_df["Strain_Type"]], ignore_index=True)
combined_df_scaled["Time"] = ["Current"]*len(df) + ["Future 10y"]*len(future_df)

# -----------------------------
# PCA for visualization
# -----------------------------
pca = PCA(n_components=2, random_state=42)
pcs = pca.fit_transform(combined_df_scaled[traits])
combined_df_scaled["PC1"] = pcs[:,0]
combined_df_scaled["PC2"] = pcs[:,1]
print("Explained variance ratio (PCA):", pca.explained_variance_ratio_)

# -----------------------------
# Plot current vs future in PCA space
# -----------------------------
plt.figure(figsize=(10,7))
sns.scatterplot(data=combined_df_scaled, x="PC1", y="PC2",
                hue="Strain_Type", style="Time", palette="tab10", s=100)
plt.title("PCA of Yeast Strains: Current vs 10-Year Projection")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(bbox_to_anchor=(1.05,1), loc="upper left")
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.animation import FuncAnimation

# -----------------------------
# Load original CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

traits = ["Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
          "Copper_Resistance","Sugar_Utilization","H2S","Protease",
          "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
          "Esterase","Gluconic_Acid_Utilization"]

# -----------------------------
# Add small noise to low-variance traits
# -----------------------------
low_variance_traits = ["H2S","Esterase","SO2_Resistance","Sugar_Utilization",
                       "Lipase","Growth_Rate","Glucosidase","Gluconic_Acid_Utilization"]
np.random.seed(42)
for t in low_variance_traits:
    if df[t].std() < 1:
        df[t] += np.random.normal(0, 0.2, size=len(df))

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
# Simulate 10-year future projection
# -----------------------------
years = 10
change_fraction = 0.05  # 5% of trait SD per year
future_df_scaled = df_scaled.copy()
for t in traits:
    std_t = df_scaled[t].std()
    yearly_change = std_t * change_fraction
    future_df_scaled[t] = df_scaled[t] + np.random.normal(0, yearly_change, size=len(df)) * years

# -----------------------------
# PCA on combined dataset
# -----------------------------
pca = PCA(n_components=2, random_state=42)
combined = pd.concat([df_scaled[traits], future_df_scaled[traits]])
pca.fit(combined)

df_scaled[["PC1","PC2"]] = pca.transform(df_scaled[traits])
future_df_scaled[["PC1","PC2"]] = pca.transform(future_df_scaled[traits])

# -----------------------------
# Animated PCA
# -----------------------------
fig, ax = plt.subplots(figsize=(10,7))
palette = {"Type strain":"tab:blue", "Environmental":"tab:green", "Clinical":"tab:red"}

# Combine for animation
years_list = np.linspace(0, 10, 50)  # 50 frames
scatter_current = ax.scatter([], [], s=100)
scatter_future = ax.scatter([], [], s=100, marker="X")

def init():
    ax.set_xlim(min(df_scaled["PC1"].min(), future_df_scaled["PC1"].min()) - 1,
                max(df_scaled["PC1"].max(), future_df_scaled["PC1"].max()) + 1)
    ax.set_ylim(min(df_scaled["PC2"].min(), future_df_scaled["PC2"].min()) - 1,
                max(df_scaled["PC2"].max(), future_df_scaled["PC2"].max()) + 1)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of Yeast Strains: 10-Year Projection")
    return scatter_current, scatter_future

def update(frame):
    # Linear interpolation between current and future
    fraction = frame / years_list[-1]
    PC1_interp = df_scaled["PC1"] + fraction * (future_df_scaled["PC1"] - df_scaled["PC1"])
    PC2_interp = df_scaled["PC2"] + fraction * (future_df_scaled["PC2"] - df_scaled["PC2"])
    ax.clear()
    for strain_type in df_scaled["Strain_Type"].unique():
        mask = df_scaled["Strain_Type"] == strain_type
        ax.scatter(PC1_interp[mask], PC2_interp[mask],
                   label=strain_type, color=palette[strain_type], s=100)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"PCA of Yeast Strains: Year {frame:.1f}")
    ax.legend()
    return ax,

ani = FuncAnimation(fig, update, frames=years_list, init_func=init, blit=False, repeat=True)
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# -----------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

traits = ["Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
          "Copper_Resistance","Sugar_Utilization","H2S","Protease",
          "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
          "Esterase","Gluconic_Acid_Utilization"]

# -----------------------------
# Add small noise to low-variance traits
# -----------------------------
low_variance_traits = ["H2S","Esterase","SO2_Resistance","Sugar_Utilization",
                       "Lipase","Growth_Rate","Glucosidase","Gluconic_Acid_Utilization"]
np.random.seed(42)
for t in low_variance_traits:
    if df[t].std() < 1:
        df[t] += np.random.normal(0, 0.2, size=len(df))

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
# Simulate 10-year future projection
# -----------------------------
years = 10
change_fraction = 0.05  # 5% of trait SD per year
future_df_scaled = df_scaled.copy()
for t in traits:
    std_t = df_scaled[t].std()
    yearly_change = std_t * change_fraction
    future_df_scaled[t] = df_scaled[t] + np.random.normal(0, yearly_change, size=len(df)) * years

# -----------------------------
# PCA on combined dataset
# -----------------------------
pca = PCA(n_components=2, random_state=42)
combined = pd.concat([df_scaled[traits], future_df_scaled[traits]])
pca.fit(combined)
df_scaled[["PC1","PC2"]] = pca.transform(df_scaled[traits])
future_df_scaled[["PC1","PC2"]] = pca.transform(future_df_scaled[traits])

# -----------------------------
# Animated PCA with moving labels
# -----------------------------
fig, ax = plt.subplots(figsize=(10,7))
palette = {"Type strain":"tab:blue", "Environmental":"tab:green", "Clinical":"tab:red"}

years_list = np.linspace(0, 10, 100)  # 100 frames for smooth animation

# Prepare scatter points
scatter_points = {strain: ax.plot([], [], 'o', color=palette[stype], markersize=8)[0]
                  for strain, stype in zip(df_scaled["Strain"], df_scaled["Strain_Type"])}
# Prepare labels
scatter_labels = {strain: ax.text(0, 0, strain, fontsize=7) for strain in df_scaled["Strain"]}

def init():
    ax.set_xlim(min(df_scaled["PC1"].min(), future_df_scaled["PC1"].min()) - 1,
                max(df_scaled["PC1"].max(), future_df_scaled["PC1"].max()) + 1)
    ax.set_ylim(min(df_scaled["PC2"].min(), future_df_scaled["PC2"].min()) - 1,
                max(df_scaled["PC2"].max(), future_df_scaled["PC2"].max()) + 1)
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of Yeast Strains: 10-Year Projection")
    return list(scatter_points.values()) + list(scatter_labels.values())

def update(frame):
    fraction = frame / years_list[-1]  # 0 -> 1
    for strain in df_scaled["Strain"]:
        # Linear interpolation between current and future
        x = df_scaled.loc[df_scaled["Strain"]==strain, "PC1"].values[0] + \
            fraction * (future_df_scaled.loc[future_df_scaled["Strain"]==strain, "PC1"].values[0] - 
                        df_scaled.loc[df_scaled["Strain"]==strain, "PC1"].values[0])
        y = df_scaled.loc[df_scaled["Strain"]==strain, "PC2"].values[0] + \
            fraction * (future_df_scaled.loc[future_df_scaled["Strain"]==strain, "PC2"].values[0] - 
                        df_scaled.loc[df_scaled["Strain"]==strain, "PC2"].values[0])
        scatter_points[strain].set_data(x, y)
        scatter_labels[strain].set_position((x+0.05, y+0.05))  # offset labels slightly
    ax.set_title(f"PCA of Yeast Strains: Year {frame:.1f}")
    return list(scatter_points.values()) + list(scatter_labels.values())

ani = FuncAnimation(fig, update, frames=years_list, init_func=init, blit=True, repeat=True)

plt.show()

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# -----------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\phen_scaled.csv"
df = pd.read_csv(csv_path)

traits = ["Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
          "Copper_Resistance","Sugar_Utilization","Protease",
          "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
          "Esterase","Gluconic_Acid_Utilization"]

# -----------------------------
# Compute correlation matrix
# -----------------------------
corr_matrix = df[traits].corr()
print("Correlation matrix:")
print(corr_matrix)

# -----------------------------
# Build phenotypic network (edges with correlation >= 0.5)
# -----------------------------
threshold = 0.5
G = nx.Graph()

# Add nodes
for trait in traits:
    G.add_node(trait)

# Add edges based on correlation threshold
for i, t1 in enumerate(traits):
    for j, t2 in enumerate(traits):
        if i < j and corr_matrix.loc[t1, t2] >= threshold:
            G.add_edge(t1, t2, weight=corr_matrix.loc[t1, t2])

# -----------------------------
# Visualize network
# -----------------------------
plt.figure(figsize=(12,10))

# Node positions
pos = nx.spring_layout(G, seed=42)

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_size=1200, node_color='skyblue')

# Draw edges with width proportional to correlation
edges = G.edges(data=True)
nx.draw_networkx_edges(G, pos, edgelist=edges, 
                       width=[d['weight']*5 for (u,v,d) in edges], alpha=0.7)

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=10, font_weight='bold')

plt.title("Phenotypic Network (Correlation >= 0.5)", fontsize=16)
plt.axis('off')
plt.show()

import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# --------------------------
# Select numeric traits
# --------------------------
numeric_cols = ["Ethanol_Tolerance", "Killer_Activity",
                "Pigment_Intensity", "Copper_Resistance"]  # adjust if needed
df_num = df[numeric_cols].copy()

# --------------------------
# Compute correlation and p-values
# --------------------------
corr_matrix = np.zeros((len(numeric_cols), len(numeric_cols)))
p_values = np.zeros((len(numeric_cols), len(numeric_cols)))

for i in range(len(numeric_cols)):
    for j in range(len(numeric_cols)):
        r, p = pearsonr(df_num.iloc[:, i], df_num.iloc[:, j])
        corr_matrix[i, j] = r
        p_values[i, j] = p

# --------------------------
# Build Network Graph (|r| ≥ 0.5 & p < 0.05)
# --------------------------
G = nx.Graph()

# Add nodes
for trait in numeric_cols:
    G.add_node(trait)

# Add edges
for i, t1 in enumerate(numeric_cols):
    for j, t2 in enumerate(numeric_cols):
        if i < j and abs(corr_matrix[i, j]) >= 0.5 and p_values[i, j] < 0.05:
            G.add_edge(t1, t2, weight=corr_matrix[i, j])

# --------------------------
# Plot the network
# --------------------------
plt.figure(figsize=(8,6))
pos = nx.spring_layout(G, seed=42)  # fixed layout

# Nodes
nx.draw_networkx_nodes(G, pos, node_color="skyblue", node_size=1500)

# Edges weighted by correlation strength
weights = [abs(G[u][v]["weight"])*3 for u,v in G.edges()]
nx.draw_networkx_edges(G, pos, width=weights)

# Labels
nx.draw_networkx_labels(G, pos, font_size=10)

plt.title("Phenotype Correlation Network (|r| ≥ 0.5 & p < 0.05)")
plt.axis("off")
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# -----------------------------
# Load CSV
# -----------------------------
csv_path = r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv"
df = pd.read_csv(csv_path)

traits = [
    "Ethanol_Tolerance","Killer_Activity","Pigment_Intensity",
    "Copper_Resistance","Sugar_Utilization","H2S","Protease",
    "Lipase","Growth_Rate","Glucosidase","SO2_Resistance",
    "Esterase","Gluconic_Acid_Utilization"
]

# -----------------------------
# Add small noise to low-variance traits
# -----------------------------
low_variance_traits = [
    "H2S","Esterase","SO2_Resistance","Sugar_Utilization",
    "Lipase","Growth_Rate","Glucosidase","Gluconic_Acid_Utilization"
]

np.random.seed(42)
for t in low_variance_traits:
    if df[t].std() < 1:
        df[t] += np.random.normal(0, 0.2, size=len(df))

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
# Simulate 10-year PAST projection
# -----------------------------
years = 10
change_fraction = 0.05  # 5% of SD per year

past_df_scaled = df_scaled.copy()
for t in traits:
    std_t = df_scaled[t].std()
    yearly_change = std_t * change_fraction
    past_df_scaled[t] = (
        df_scaled[t]
        - np.random.normal(0, yearly_change, size=len(df)) * years
    )

# -----------------------------
# PCA on combined (past + present)
# -----------------------------
pca = PCA(n_components=2, random_state=42)
combined = pd.concat([past_df_scaled[traits], df_scaled[traits]])
pca.fit(combined)

past_df_scaled[["PC1","PC2"]] = pca.transform(past_df_scaled[traits])
df_scaled[["PC1","PC2"]] = pca.transform(df_scaled[traits])

# -----------------------------
# Animated PCA (Past → Present)
# -----------------------------
fig, ax = plt.subplots(figsize=(10,7))

palette = {
    "Type strain": "tab:blue",
    "Environmental": "tab:green",
    "Clinical": "tab:red"
}

frames = np.linspace(0, 10, 100)  # -10 → 0 years

scatter_points = {
    strain: ax.plot([], [], 'o',
                    color=palette[stype],
                    markersize=8)[0]
    for strain, stype in zip(df_scaled["Strain"], df_scaled["Strain_Type"])
}

scatter_labels = {
    strain: ax.text(0, 0, strain, fontsize=7)
    for strain in df_scaled["Strain"]
}

def init():
    ax.set_xlim(
        min(past_df_scaled["PC1"].min(), df_scaled["PC1"].min()) - 1,
        max(past_df_scaled["PC1"].max(), df_scaled["PC1"].max()) + 1
    )
    ax.set_ylim(
        min(past_df_scaled["PC2"].min(), df_scaled["PC2"].min()) - 1,
        max(past_df_scaled["PC2"].max(), df_scaled["PC2"].max()) + 1
    )
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PCA of Yeast Strains: 10 Years Before")
    return list(scatter_points.values()) + list(scatter_labels.values())

def update(frame):
    fraction = frame / frames[-1]  # 0 → 1
    year = -10 + frame

    for strain in df_scaled["Strain"]:
        x0 = past_df_scaled.loc[past_df_scaled["Strain"]==strain, "PC1"].values[0]
        y0 = past_df_scaled.loc[past_df_scaled["Strain"]==strain, "PC2"].values[0]
        x1 = df_scaled.loc[df_scaled["Strain"]==strain, "PC1"].values[0]
        y1 = df_scaled.loc[df_scaled["Strain"]==strain, "PC2"].values[0]

        x = x0 + fraction * (x1 - x0)
        y = y0 + fraction * (y1 - y0)

        scatter_points[strain].set_data(x, y)
        scatter_labels[strain].set_position((x+0.05, y+0.05))

    ax.set_title(f"PCA of Yeast Strains: Year {year:.1f}")
    return list(scatter_points.values()) + list(scatter_labels.values())

ani = FuncAnimation(
    fig,
    update,
    frames=frames,
    init_func=init,
    blit=True,
    repeat=True
)

plt.show()

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
