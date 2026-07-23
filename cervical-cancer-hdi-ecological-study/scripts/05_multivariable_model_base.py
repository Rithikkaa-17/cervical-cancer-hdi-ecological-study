import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv")
df = df.dropna(subset=["Screening_program_bin"]).copy()
print("Analytic n =", len(df))

df["logGNI"] = np.log(df["GNI_per_capita"])
df["logASIR"] = np.log(df["ASIR"])
df["logASMR"] = np.log(df["ASMR"])

# Bivariate correlations
print("\n=== BIVARIATE SPEARMAN CORRELATIONS ===")
for var in ["HDI","HPV_vax_coverage","GNI_per_capita","Screening_program_bin"]:
    for outcome in ["ASIR","ASMR"]:
        rho,p = stats.spearmanr(df[var], df[outcome])
        print(f"{var} vs {outcome}: rho={rho:.3f} p={p:.2e}")

# VIF check
X_vif = df[["HDI","HPV_vax_coverage","logGNI","Screening_program_bin"]].dropna()
X_vif_c = sm.add_constant(X_vif)
print("\n=== VIF ===")
for i, col in enumerate(X_vif_c.columns):
    if col == "const": continue
    vif = variance_inflation_factor(X_vif_c.values, i)
    print(f"{col}: VIF={vif:.2f}")

# Multivariable regression
predictors = ["HDI","HPV_vax_coverage","logGNI","Screening_program_bin"]
X = sm.add_constant(df[predictors])

print("\n=== OLS: log(ASIR) ~ HDI + HPV_vax + logGNI + Screening ===")
m1 = sm.OLS(df["logASIR"], X).fit()
print(m1.summary().tables[1])
print("R2 =", m1.rsquared)

print("\n=== OLS: log(ASMR) ~ HDI + HPV_vax + logGNI + Screening ===")
m2 = sm.OLS(df["logASMR"], X).fit()
print(m2.summary().tables[1])
print("R2 =", m2.rsquared)

df.to_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", index=False)

print("\n" + "="*70)
print("VIF was high for HDI & logGNI (collinear, since GNI feeds into HDI).")
print("Re-running without GNI to get stable estimates for HDI, HPV vax, Screening.")
print("="*70)

predictors2 = ["HDI","HPV_vax_coverage","Screening_program_bin"]
X2 = sm.add_constant(df[predictors2])

X2_vif = sm.add_constant(df[predictors2])
print("\n=== VIF (reduced model) ===")
for i, col in enumerate(X2_vif.columns):
    if col == "const": continue
    vif = variance_inflation_factor(X2_vif.values, i)
    print(f"{col}: VIF={vif:.2f}")

print("\n=== OLS: log(ASIR) ~ HDI + HPV_vax + Screening ===")
m3 = sm.OLS(df["logASIR"], X2).fit()
print(m3.summary().tables[1])
print("R2 =", m3.rsquared)

print("\n=== OLS: log(ASMR) ~ HDI + HPV_vax + Screening ===")
m4 = sm.OLS(df["logASMR"], X2).fit()
print(m4.summary().tables[1])
print("R2 =", m4.rsquared)
