import pandas as pd

# Load merged cervical cancer data
cc = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_globocan2024_merged.csv", keep_default_na=False, na_values=[""])
print("Cervical cancer data:", cc.shape)

# Load HDI data
hdi = pd.read_csv("/mnt/user-data/uploads/human-development-index.csv", keep_default_na=False, na_values=[""])
hdi.columns = ["Country", "Alpha-3 code", "Year", "HDI", "OWID_region"]

# Drop aggregate rows (no ISO-3 code = regions/income groups, not countries)
hdi_countries = hdi[hdi["Alpha-3 code"].notna() & (hdi["Alpha-3 code"] != "")].copy()
print("After dropping aggregates:", hdi_countries.shape, "| unique countries:", hdi_countries["Alpha-3 code"].nunique())

# Keep latest year per country
hdi_latest = hdi_countries.sort_values("Year").groupby("Alpha-3 code", as_index=False).last()
hdi_latest = hdi_latest[["Alpha-3 code", "Country", "Year", "HDI"]].rename(columns={"Year": "HDI_year", "Country": "HDI_country_name"})
print("Latest-year HDI:", hdi_latest.shape)
print(hdi_latest["HDI_year"].value_counts())

# Merge with cervical cancer data
merged = cc.merge(hdi_latest, on="Alpha-3 code", how="left")
print("\nFinal merged shape:", merged.shape)
print("Countries missing HDI:", merged["HDI"].isna().sum())
print(merged.loc[merged["HDI"].isna(), "Country"].tolist())

merged.to_csv("/mnt/user-data/outputs/cervical_cancer_with_HDI.csv", index=False)
