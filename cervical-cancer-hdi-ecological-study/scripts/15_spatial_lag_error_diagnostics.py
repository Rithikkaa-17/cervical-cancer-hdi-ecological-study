import pandas as pd
import numpy as np
from libpysal.weights import KNN
from spreg import ML_Lag, ML_Error, OLS

df = pd.read_csv("cervical_cancer_spatial.csv")
df["logASIR"] = np.log(df["ASIR"])

predictors = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female","HIV_prev"]
y = df[["logASIR"]].values
X = df[predictors].values

coords = df[["longitude","latitude"]].values
w = KNN.from_array(coords, k=5)
w.transform = "r"

# --- Standard OLS (for comparison) ---
ols = OLS(y, X, w=w, name_y="logASIR", name_x=predictors, spat_diag=True, moran=True)
print("="*70)
print("STANDARD OLS (spatial diagnostics)")
print("="*70)
print(f"R2: {ols.r2:.4f}")
for name, coef, se, t, p in zip(["const"]+predictors, ols.betas.flatten(), np.sqrt(np.diag(ols.vm)), ols.t_stat, ols.t_stat):
    pass
print(ols.summary[ols.summary.find("REGRESSION"):ols.summary.find("REGRESSION")+3000])
