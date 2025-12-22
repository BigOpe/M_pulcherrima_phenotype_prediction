import pandas as pd
import pingouin as pg

# Load your CSV
df = pd.read_csv(r"C:\Users\Giorgia\dna\dna2vec-legacy\Scripts\phen.csv")

# Set strain as index (recommended for clarity)
df = df.set_index("Strain")

# Keep only numeric phenotypic traits
df_numeric = df.drop(columns=["Notes"])

# Double-check all columns are numeric
print(df_numeric.dtypes)
df_numeric = df_numeric.apply(pd.to_numeric, errors="coerce")
df_numeric = df_numeric.fillna(df_numeric.mean())
alpha, ci = pg.cronbach_alpha(df_numeric)

print(f"Cronbach’s alpha: {alpha:.3f}")
print(f"95% CI: {ci}")
