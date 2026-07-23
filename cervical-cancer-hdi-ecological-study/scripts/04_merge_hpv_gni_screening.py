import pandas as pd
import numpy as np

def latest_per_country(path, value_col, out_name, drop_na_value=True):
    df = pd.read_csv(f"/mnt/user-data/uploads/{path}", keep_default_na=False, na_values=[""])
    df = df[df["code"].notna() & (df["code"] != "")].copy()
    df["year"] = pd.to_numeric(df["year"], errors="coerce")
    if drop_na_value:
        df[value_col] = pd.to_numeric(df[value_col], errors="coerce") if df[value_col].dtype == object else df[value_col]
    df = df.sort_values("year")
    latest = df.groupby("code", as_index=False).last()
    return latest[["code", "year", value_col]].rename(columns={"year": f"{out_name}_year", value_col: out_name})

# Base: cervical cancer + HDI (already merged)
base = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_with_HDI.csv", keep_default_na=False, na_values=[""])
base = base.rename(columns={"Alpha-3 code": "code"})

# HPV vaccination coverage
hpv = latest_per_country("coverage-of-the-human-papillomavirus-vaccine.csv",
                          "coverage__antigen_description_hpv_vaccination_coverage_by_age_15__last_dose__females",
                          "HPV_vax_coverage")
print("HPV rows:", hpv.shape, "| non-null:", hpv["HPV_vax_coverage"].notna().sum())

# GNI per capita
gni = latest_per_country("gross-national-income-per-capita-undp.csv", "gni_pc__sex_total", "GNI_per_capita")
print("GNI rows:", gni.shape, "| non-null:", gni["GNI_per_capita"].notna().sum())

# Screening program (categorical Yes/No)
scr = pd.read_csv("/mnt/user-data/uploads/countries-with-national-cervical-cancer-screening-program.csv", keep_default_na=False, na_values=[""])
scr = scr[scr["code"].notna() & (scr["code"] != "")].copy()
scr["year"] = pd.to_numeric(scr["year"], errors="coerce")
scr = scr.sort_values("year").groupby("code", as_index=False).last()
scr = scr[["code", "year", "existence_of_national_screening_program_for_cervical_cancer"]].rename(
    columns={"year": "Screening_year", "existence_of_national_screening_program_for_cervical_cancer": "Screening_program"})
scr["Screening_program_bin"] = (scr["Screening_program"] == "Yes").astype(int)
print("Screening rows:", scr.shape, "| Yes count:", scr["Screening_program_bin"].sum())

# Merge all
merged = base.merge(hpv[["code","HPV_vax_coverage"]], on="code", how="left")
merged = merged.merge(gni[["code","GNI_per_capita"]], on="code", how="left")
merged = merged.merge(scr[["code","Screening_program","Screening_program_bin"]], on="code", how="left")

# HPV coverage: missing -> assume 0 (no vaccination program / not reported = negligible coverage)
merged["HPV_vax_coverage"] = merged["HPV_vax_coverage"].fillna(0)

print("\nFinal shape:", merged.shape)
print("Missing GNI:", merged["GNI_per_capita"].isna().sum())
print("Missing Screening:", merged["Screening_program_bin"].isna().sum())
print(merged.loc[merged["GNI_per_capita"].isna(), "Country"].tolist())
print(merged.loc[merged["Screening_program_bin"].isna(), "Country"].tolist())

merged.to_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", index=False)
