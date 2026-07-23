import pandas as pd
import numpy as np
from scipy import stats

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_with_HDI.csv")

print("="*60)
print("DESCRIPTIVE STATS (n=%d countries)" % len(df))
print("="*60)
for col in ["ASIR","ASMR","MIR","HDI"]:
    print(f"{col}: median={df[col].median():.3f}, IQR=[{df[col].quantile(.25):.3f}-{df[col].quantile(.75):.3f}], range=[{df[col].min():.3f}-{df[col].max():.3f}]")

print("\n" + "="*60)
print("SPEARMAN CORRELATIONS (HDI vs outcomes)")
print("="*60)
for outcome in ["ASIR","ASMR","MIR"]:
    rho, p = stats.spearmanr(df["HDI"], df[outcome])
    print(f"HDI vs {outcome}: rho={rho:.3f}, p={p:.2e}")

print("\n" + "="*60)
print("PEARSON CORRELATIONS (for comparison to prior literature)")
print("="*60)
for outcome in ["ASIR","ASMR","MIR"]:
    r, p = stats.pearsonr(df["HDI"], df[outcome])
    print(f"HDI vs {outcome}: r={r:.3f}, p={p:.2e}")

# Simple linear regression: log(ASIR) ~ HDI  (rates are right-skewed)
import statsmodels.api as sm

df["logASIR"] = np.log(df["ASIR"])
df["logASMR"] = np.log(df["ASMR"])

X = sm.add_constant(df["HDI"])
for outcome, label in [("logASIR","log(ASIR)"), ("logASMR","log(ASMR)")]:
    model = sm.OLS(df[outcome], X).fit()
    print(f"\n--- OLS: {label} ~ HDI ---")
    print(f"HDI coef = {model.params['HDI']:.3f} (p={model.pvalues['HDI']:.2e}), R^2={model.rsquared:.3f}")

# HDI tiers (UNDP standard cutoffs)
def hdi_tier(h):
    if h >= 0.8: return "Very High"
    elif h >= 0.7: return "High"
    elif h >= 0.55: return "Medium"
    else: return "Low"

df["HDI_tier"] = df["HDI"].apply(hdi_tier)
tier_order = ["Low","Medium","High","Very High"]
print("\n" + "="*60)
print("BY HDI TIER")
print("="*60)
summary = df.groupby("HDI_tier")[["ASIR","ASMR","MIR"]].median().reindex(tier_order)
counts = df["HDI_tier"].value_counts().reindex(tier_order)
summary["n_countries"] = counts
print(summary)

print("\nRatio Low:VeryHigh ASIR =", round(summary.loc["Low","ASIR"]/summary.loc["Very High","ASIR"],2))
print("Ratio Low:VeryHigh ASMR =", round(summary.loc["Low","ASMR"]/summary.loc["Very High","ASMR"],2))

df.to_csv("/mnt/user-data/outputs/cervical_cancer_with_HDI.csv", index=False)
