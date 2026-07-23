import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv")
df = df.dropna(subset=["Screening_program_bin","Smoking_prev_female"]).copy()
print("Analytic n =", len(df))

df["logASIR"] = np.log(df["ASIR"])
df["logASMR"] = np.log(df["ASMR"])

print("\n=== BIVARIATE: Smoking vs outcomes ===")
for outcome in ["ASIR","ASMR"]:
    rho,p = stats.spearmanr(df["Smoking_prev_female"], df[outcome])
    print(f"Smoking vs {outcome}: rho={rho:.3f} p={p:.2e}")

predictors = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female"]
X = sm.add_constant(df[predictors])

print("\n=== VIF ===")
for i, col in enumerate(X.columns):
    if col == "const": continue
    vif = variance_inflation_factor(X.values, i)
    print(f"{col}: VIF={vif:.2f}")

print("\n=== OLS: log(ASIR) ~ HDI + HPV_vax + Screening + Smoking ===")
m1 = sm.OLS(df["logASIR"], X).fit()
print(m1.summary().tables[1])
print("R2 =", m1.rsquared)

print("\n=== OLS: log(ASMR) ~ HDI + HPV_vax + Screening + Smoking ===")
m2 = sm.OLS(df["logASMR"], X).fit()
print(m2.summary().tables[1])
print("R2 =", m2.rsquared)

df.to_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv", index=False)
