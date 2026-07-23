import pandas as pd

inc = pd.read_csv("/mnt/user-data/uploads/dataset-inc-females-in-2024-cervix-uteri.csv", keep_default_na=False, na_values=[""])
mort = pd.read_csv("/mnt/user-data/uploads/dataset-mort-females-in-2024-cervix-uteri.csv", keep_default_na=False, na_values=[""])

print("Incidence shape:", inc.shape)
print("Mortality shape:", mort.shape)
print(inc.columns.tolist())

inc_s = inc[["Alpha-3 code", "Label", "Number", "ASR (World)", "Crude rate", "Cumulative risk"]].rename(
    columns={"Label": "Country", "Number": "Cases", "ASR (World)": "ASIR", "Crude rate": "Crude_incidence_rate", "Cumulative risk": "Cum_risk_incidence"}
)
mort_s = mort[["Alpha-3 code", "Number", "ASR (World)", "Crude rate", "Cumulative risk"]].rename(
    columns={"Number": "Deaths", "ASR (World)": "ASMR", "Crude rate": "Crude_mortality_rate", "Cumulative risk": "Cum_risk_mortality"}
)

merged = inc_s.merge(mort_s, on="Alpha-3 code", how="outer")
merged["ASIR"] = pd.to_numeric(merged["ASIR"], errors="coerce")
merged["ASMR"] = pd.to_numeric(merged["ASMR"], errors="coerce")
merged["MIR"] = (merged["ASMR"] / merged["ASIR"]).round(3)

merged = merged.sort_values("ASIR", ascending=False)
print(merged.shape)
print(merged.head(10))
print("Missing ASIR:", merged["ASIR"].isna().sum(), "Missing ASMR:", merged["ASMR"].isna().sum())

merged.to_csv("/mnt/user-data/outputs/cervical_cancer_globocan2024_merged.csv", index=False)
