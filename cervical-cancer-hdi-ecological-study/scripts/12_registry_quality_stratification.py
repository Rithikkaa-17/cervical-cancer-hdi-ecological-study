import pandas as pd
import numpy as np
from scipy import stats

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv")
df_reg = df.dropna(subset=["Death_reg_pct"]).copy()
print(f"n with registration data: {len(df_reg)}")
print(f"Death registration completeness: median={df_reg['Death_reg_pct'].median():.1f}%, "
      f"IQR=[{df_reg['Death_reg_pct'].quantile(.25):.1f}-{df_reg['Death_reg_pct'].quantile(.75):.1f}]")

# Split at 90% completeness -- a standard "high-quality vital statistics" threshold used in demography literature
threshold = 90
df_reg["registry_quality"] = np.where(df_reg["Death_reg_pct"] >= threshold, "High (>=90%)", "Lower (<90%)")
print(df_reg["registry_quality"].value_counts())

print("\n" + "="*65)
print("HDI-ASIR CORRELATION, STRATIFIED BY REGISTRATION QUALITY")
print("="*65)
for grp in ["High (>=90%)", "Lower (<90%)"]:
    sub = df_reg[df_reg["registry_quality"]==grp]
    rho_i, p_i = stats.spearmanr(sub["HDI"], sub["ASIR"])
    rho_m, p_m = stats.spearmanr(sub["HDI"], sub["ASMR"])
    print(f"\n{grp} (n={len(sub)}):")
    print(f"  HDI vs ASIR: rho={rho_i:.3f}, p={p_i:.2e}")
    print(f"  HDI vs ASMR: rho={rho_m:.3f}, p={p_m:.2e}")

# Also check: is HDI itself correlated with registration completeness? (expected, but quantify)
rho_check, p_check = stats.spearmanr(df_reg["HDI"], df_reg["Death_reg_pct"])
print(f"\nHDI vs Death_reg_pct: rho={rho_check:.3f}, p={p_check:.2e} (expected strong positive - registry quality tracks development)")

# Tier-level medians within each stratum
print("\n" + "="*65)
print("MEDIAN ASIR BY HDI TIER, WITHIN EACH REGISTRATION-QUALITY STRATUM")
print("="*65)
tier_order = ["Low","Medium","High","Very High"]
for grp in ["High (>=90%)", "Lower (<90%)"]:
    sub = df_reg[df_reg["registry_quality"]==grp]
    print(f"\n{grp}:")
    summary = sub.groupby("HDI_tier")["ASIR"].agg(["median","count"]).reindex(tier_order)
    print(summary)

df_reg.to_csv("/mnt/user-data/outputs/cervical_cancer_registry_quality.csv", index=False)
