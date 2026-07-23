import pandas as pd
import numpy as np

# --- 2024 cycle: our existing full dataset ---
df2024 = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", keep_default_na=False, na_values=[""])
df2024 = df2024.rename(columns={"code":"Code"})[["Code","Country","ASIR","ASMR","HDI"]].copy()
df2024["Year"] = 2024
df2024 = df2024.rename(columns={"ASIR":"ASIR_val","ASMR":"ASMR_val","HDI":"HDI_val"})

# --- 2022 cycle: GLOBOCAN 2022 incidence (via OWID) + time-matched HDI (2021, closest available year) ---
df2022_inc = pd.read_csv("globocan2022.csv", keep_default_na=False, na_values=[""])
df2022_inc = df2022_inc.rename(columns={"Entity":"Country"})

hdi_full = pd.read_csv("human-development-index.csv", keep_default_na=False, na_values=[""])
hdi_full.columns = ["Country","Code","Year","HDI","Region"]
hdi_full = hdi_full[hdi_full["Code"].notna() & (hdi_full["Code"]!="")].copy()
hdi_full["Year"] = pd.to_numeric(hdi_full["Year"], errors="coerce")
hdi_2021 = hdi_full[hdi_full["Year"]==2021][["Code","HDI"]].rename(columns={"HDI":"HDI_val"})

df2022 = df2022_inc.merge(hdi_2021, on="Code", how="left")
df2022 = df2022.rename(columns={"ASIR_2022":"ASIR_val"})
df2022["ASMR_val"] = np.nan  # OWID GLOBOCAN-2022 mirror only has incidence, not mortality
df2022["Year"] = 2022
df2022 = df2022[["Code","Country","ASIR_val","ASMR_val","HDI_val","Year"]]

print("2022 cycle:", df2022.shape, "| missing HDI:", df2022["HDI_val"].isna().sum())
print("2024 cycle:", df2024.shape, "| missing HDI:", df2024["HDI_val"].isna().sum())

# --- Combine into panel (long format) ---
panel = pd.concat([df2022, df2024[["Code","Country","ASIR_val","ASMR_val","HDI_val","Year"]]], ignore_index=True)
panel = panel.dropna(subset=["HDI_val","ASIR_val"])
print("\nCombined panel:", panel.shape)
print(panel["Year"].value_counts())

panel.to_csv("cervical_cancer_panel_2022_2024.csv", index=False)
