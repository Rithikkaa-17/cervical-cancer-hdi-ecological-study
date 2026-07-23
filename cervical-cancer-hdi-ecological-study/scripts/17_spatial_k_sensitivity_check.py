import pandas as pd
import numpy as np
from libpysal.weights import KNN
from esda.moran import Moran
from spreg import ML_Lag
import statsmodels.api as sm

df = pd.read_csv("cervical_cancer_spatial.csv")
df["logASIR"] = np.log(df["ASIR"])
predictors = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female","HIV_prev"]
y = df[["logASIR"]].values
X = df[predictors].values
coords = df[["longitude","latitude"]].values

print(f"{'k':<5}{'Moran I (resid)':<18}{'p-value':<12}{'HDI direct (lag)':<20}{'HDI p-value':<12}{'Spatial rho':<14}{'rho p-value'}")
for k in [3, 5, 8, 10]:
    w = KNN.from_array(coords, k=k)
    w.transform = "r"

    # OLS residuals for Moran's I
    Xc = sm.add_constant(df[predictors])
    ols = sm.OLS(df["logASIR"], Xc).fit()
    mi = Moran(ols.resid.values, w)

    # Spatial lag model
    lag = ML_Lag(y, X, w=w, name_y="logASIR", name_x=predictors)
    hdi_idx = predictors.index("HDI")
    hdi_coef = lag.betas[hdi_idx+1][0]  # +1 for constant
    hdi_se = np.sqrt(lag.vm[hdi_idx+1, hdi_idx+1])
    hdi_z = hdi_coef / hdi_se
    from scipy import stats as st
    hdi_p = 2*(1-st.norm.cdf(abs(hdi_z)))
    rho = lag.betas[-1][0]
    rho_se = np.sqrt(lag.vm[-1,-1])
    rho_p = 2*(1-st.norm.cdf(abs(rho/rho_se)))

    print(f"{k:<5}{mi.I:<18.3f}{mi.p_sim:<12.4f}{hdi_coef:<20.3f}{hdi_p:<12.2e}{rho:<14.3f}{rho_p:.2e}")
