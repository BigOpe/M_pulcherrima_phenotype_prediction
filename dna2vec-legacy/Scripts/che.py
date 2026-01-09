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
