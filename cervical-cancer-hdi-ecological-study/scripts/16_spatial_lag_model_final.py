import pandas as pd
import numpy as np
from libpysal.weights import KNN
from spreg import ML_Lag

df = pd.read_csv("cervical_cancer_spatial.csv")
df["logASIR"] = np.log(df["ASIR"])

predictors = ["HDI","HPV_vax_coverage","Screening_program_bin","Smoking_prev_female","HIV_prev"]
y = df[["logASIR"]].values
X = df[predictors].values

coords = df[["longitude","latitude"]].values
w = KNN.from_array(coords, k=5)
w.transform = "r"

lag_model = ML_Lag(y, X, w=w, name_y="logASIR", name_x=predictors)
print(lag_model.summary[lag_model.summary.find("REGRESSION"):lag_model.summary.find("REGRESSION")+2200])
