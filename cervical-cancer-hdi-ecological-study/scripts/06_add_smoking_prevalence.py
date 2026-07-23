import pandas as pd
import numpy as np

base = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", keep_default_na=False, na_values=[""])

smk = pd.read_csv("/mnt/user-data/uploads/share-of-women-who-are-smoking.csv", keep_default_na=False, na_values=[""])
smk = smk[smk["code"].notna() & (smk["code"] != "")].copy()
smk["year"] = pd.to_numeric(smk["year"], errors="coerce")
smk = smk.sort_values("year").groupby("code", as_index=False).last()
smk = smk[["code","year","tobacco_use_pct_age_std__sex_female"]].rename(
    columns={"year":"Smoking_year","tobacco_use_pct_age_std__sex_female":"Smoking_prev_female"})
print("Smoking data:", smk.shape, "| non-null:", smk["Smoking_prev_female"].notna().sum())

merged = base.merge(smk[["code","Smoking_prev_female"]], on="code", how="left")
print("Merged shape:", merged.shape)
print("Missing smoking:", merged["Smoking_prev_female"].isna().sum())
print(merged.loc[merged["Smoking_prev_female"].isna(), "Country"].tolist())

merged.to_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", index=False)
