import pandas as pd
import numpy as np

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", keep_default_na=False, na_values=[""])

# --- HIV prevalence: latest available year per country ---
hiv = pd.read_csv("/mnt/user-data/uploads/share-of-the-population-infected-with-hiv.csv", keep_default_na=False, na_values=[""])
hiv.columns = ["Country","code","Year","HIV_prev"]
hiv = hiv[hiv["code"].notna() & (hiv["code"]!="")].copy()
hiv["Year"] = pd.to_numeric(hiv["Year"], errors="coerce")
hiv["HIV_prev"] = pd.to_numeric(hiv["HIV_prev"], errors="coerce")
hiv_latest = hiv.sort_values("Year").groupby("code", as_index=False).last()[["code","HIV_prev","Year"]].rename(columns={"Year":"HIV_year"})
print("HIV data:", hiv_latest.shape, "| non-null:", hiv_latest["HIV_prev"].notna().sum())

# --- Death registration completeness: average over available years (2015-2019 window) as stable estimate ---
reg = pd.read_csv("/mnt/user-data/uploads/share-of-deaths-registered.csv", keep_default_na=False, na_values=[""])
reg.columns = ["Country","code","Year","Death_reg_pct","Region"]
reg = reg[reg["code"].notna() & (reg["code"]!="")].copy()
reg["Death_reg_pct"] = pd.to_numeric(reg["Death_reg_pct"], errors="coerce")
reg_avg = reg.groupby("code", as_index=False)["Death_reg_pct"].mean()
print("Registration completeness data:", reg_avg.shape, "| non-null:", reg_avg["Death_reg_pct"].notna().sum())

merged = df.merge(hiv_latest[["code","HIV_prev"]], on="code", how="left")
merged = merged.merge(reg_avg, on="code", how="left")

print("\nFinal shape:", merged.shape)
print("Missing HIV:", merged["HIV_prev"].isna().sum())
print("Missing Death_reg_pct:", merged["Death_reg_pct"].isna().sum())

merged.to_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", index=False)
