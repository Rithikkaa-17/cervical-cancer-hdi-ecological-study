import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv")

# ============================================================
# PART 1: Add HIV to the multivariable model
# ============================================================
df_hiv = df.dropna(subset=["Screening_program_bin","Smoking_prev_female","HIV_prev"]).copy()
print(f"HIV-inclusive analytic sample: n={len(df_hiv)}")

df_hiv["logASIR"] = np.log(df_hiv["ASIR"])
df_hiv["logASMR"] = np.log(df_hiv["ASMR"])
df_hiv["logHIV"] = np.log(df_hiv["HIV_prev"] + 0.01)  # small offset for zeros

print("\nBivariate: HIV vs outcomes")
for outcome in ["ASIR","ASMR"]:
    rho, p = stats.spearmanr(df_hiv["HIV_prev"], df_hiv[outcome])
    print(f"HIV vs {outcome}: rho={rho:.3f}, p={p:.2e}")

predictors_hiv = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female","HIV_prev"]
X = sm.add_constant(df_hiv[predictors_hiv])

print("\nVIF:")
for i, col in enumerate(X.columns):
    if col=="const": continue
    print(f"  {col}: {variance_inflation_factor(X.values, i):.2f}")

print("\n=== log(ASIR) ~ HDI + HPV_vax + Screening + Smoking + HIV ===")
m1 = sm.OLS(df_hiv["logASIR"], X).fit()
print(m1.summary().tables[1])
print("R2 =", m1.rsquared)

print("\n=== log(ASMR) ~ HDI + HPV_vax + Screening + Smoking + HIV ===")
m2 = sm.OLS(df_hiv["logASMR"], X).fit()
print(m2.summary().tables[1])
print("R2 =", m2.rsquared)

# Check: does HDI coefficient shrink once HIV is added? (compare to prior model without HIV)
predictors_no_hiv = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female"]
X0 = sm.add_constant(df_hiv[predictors_no_hiv])
m0 = sm.OLS(df_hiv["logASIR"], X0).fit()
print(f"\nHDI coefficient WITHOUT HIV (same n={len(df_hiv)}): {m0.params['HDI']:.3f}")
print(f"HDI coefficient WITH HIV: {m1.params['HDI']:.3f}")
print(f"% attenuation: {100*(1 - m1.params['HDI']/m0.params['HDI']):.1f}%")

df_hiv.to_csv("/mnt/user-data/outputs/cervical_cancer_hiv_sample.csv", index=False)
