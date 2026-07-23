import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy import stats

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_hiv_sample.csv")
df["logASMR"] = np.log(df["ASMR"])
n = len(df)
print(f"n = {n}")

# Baron & Kenny / product-of-coefficients mediation: does HIV mediate HDI's effect on ASMR?
# Path a: HDI -> HIV
Xa = sm.add_constant(df["HDI"])
model_a = sm.OLS(df["HIV_prev"], Xa).fit()
a = model_a.params["HDI"]
se_a = model_a.bse["HDI"]
print(f"\nPath a (HDI -> HIV): a={a:.3f}, se={se_a:.3f}, p={model_a.pvalues['HDI']:.2e}")

# Path b: HIV -> ASMR, controlling for HDI (and other covariates for full adjustment)
Xb = sm.add_constant(df[["HDI","HIV_prev","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female"]])
model_b = sm.OLS(df["logASMR"], Xb).fit()
b = model_b.params["HIV_prev"]
se_b = model_b.bse["HIV_prev"]
c_prime = model_b.params["HDI"]  # direct effect of HDI, controlling for HIV (and others)
print(f"Path b (HIV -> logASMR | HDI): b={b:.4f}, se={se_b:.4f}, p={model_b.pvalues['HIV_prev']:.2e}")
print(f"Direct effect of HDI on logASMR (c', controlling for HIV+covariates): {c_prime:.3f}, p={model_b.pvalues['HDI']:.2e}")

# Total effect c: HDI -> ASMR, NOT controlling for HIV (bivariate-ish, but adjust for other covariates minus HIV for apples-to-apples)
Xc = sm.add_constant(df[["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female"]])
model_c = sm.OLS(df["logASMR"], Xc).fit()
c = model_c.params["HDI"]
print(f"Total effect of HDI on logASMR (c, without HIV in model): {c:.3f}, p={model_c.pvalues['HDI']:.2e}")

# Indirect effect = a*b
indirect = a * b
print(f"\nIndirect effect (a*b): {indirect:.4f}")
print(f"Direct effect (c'): {c_prime:.4f}")
print(f"Total effect (c = a*b + c'): {indirect + c_prime:.4f} (should ~= {c:.4f})")
print(f"Proportion mediated: {100*indirect/c:.1f}%")

# Sobel test for significance of indirect effect
sobel_se = np.sqrt(b**2 * se_a**2 + a**2 * se_b**2)
sobel_z = indirect / sobel_se
sobel_p = 2 * (1 - stats.norm.cdf(abs(sobel_z)))
print(f"\nSobel test: z={sobel_z:.3f}, p={sobel_p:.4f}")

# Bootstrap confidence interval for indirect effect (bias-corrected percentile, 5000 reps)
np.random.seed(42)
n_boot = 5000
boot_indirect = []
for _ in range(n_boot):
    idx = np.random.choice(n, n, replace=True)
    dboot = df.iloc[idx]
    Xa_b = sm.add_constant(dboot["HDI"])
    ma = sm.OLS(dboot["HIV_prev"], Xa_b).fit()
    Xb_b = sm.add_constant(dboot[["HDI","HIV_prev","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female"]])
    mb = sm.OLS(dboot["logASMR"], Xb_b).fit()
    boot_indirect.append(ma.params["HDI"] * mb.params["HIV_prev"])
boot_indirect = np.array(boot_indirect)
ci_lower, ci_upper = np.percentile(boot_indirect, [2.5, 97.5])
print(f"Bootstrap 95% CI for indirect effect: [{ci_lower:.4f}, {ci_upper:.4f}]")
print(f"CI excludes zero: {ci_lower > 0 or ci_upper < 0}")
