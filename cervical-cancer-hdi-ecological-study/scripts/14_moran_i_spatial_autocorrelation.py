import pandas as pd
import numpy as np
import json
import statsmodels.api as sm
from libpysal.weights import KNN
from esda.moran import Moran

# --- Load centroids ---
with open("country_centroids.json") as f:
    centroids = json.load(f)
coords_df = pd.DataFrame(centroids).rename(columns={"alpha3":"code"})
print("Centroids loaded:", coords_df.shape)

# --- Load our HIV-inclusive analytic dataset (n=131) ---
df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_hiv_sample.csv")
df = df.merge(coords_df[["code","latitude","longitude"]], on="code", how="left")
print("Merged with coords:", df.shape, "| missing coords:", df["latitude"].isna().sum())
print(df.loc[df["latitude"].isna(), "Country"].tolist())

df = df.dropna(subset=["latitude","longitude"]).reset_index(drop=True)
print("Final n for spatial analysis:", len(df))

# --- Fit the core multivariable model (same as manuscript Section 3.6) ---
df["logASIR"] = np.log(df["ASIR"])
predictors = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female","HIV_prev"]
X = sm.add_constant(df[predictors])
model = sm.OLS(df["logASIR"], X).fit()
df["residuals"] = model.resid
print("\nModel R2:", model.rsquared)

# --- Build k-nearest-neighbor spatial weights (k=5, standard choice for country-level analysis) ---
coords = df[["longitude","latitude"]].values
w = KNN.from_array(coords, k=5)
w.transform = "r"  # row-standardize

# --- Moran's I on raw ASIR (is incidence itself spatially clustered?) ---
moran_raw = Moran(df["ASIR"].values, w)
print(f"\nMoran's I on raw ASIR: I={moran_raw.I:.4f}, p={moran_raw.p_sim:.4f} (999 permutations)")

# --- Moran's I on model residuals (is what the model CAN'T explain still spatially clustered?) ---
moran_resid = Moran(df["residuals"].values, w)
print(f"Moran's I on model residuals: I={moran_resid.I:.4f}, p={moran_resid.p_sim:.4f} (999 permutations)")

# --- Moran's I on HDI itself, for context ---
moran_hdi = Moran(df["HDI"].values, w)
print(f"Moran's I on HDI (for context): I={moran_hdi.I:.4f}, p={moran_hdi.p_sim:.4f}")

df.to_csv("cervical_cancer_spatial.csv", index=False)
