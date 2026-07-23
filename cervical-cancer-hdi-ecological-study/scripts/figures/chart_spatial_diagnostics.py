import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from libpysal.weights import KNN
from esda.moran import Moran

df = pd.read_csv("cervical_cancer_spatial.csv")
coords = df[["longitude","latitude"]].values
w = KNN.from_array(coords, k=5)
w.transform = "r"
moran_resid = Moran(df["residuals"].values, w)

PALETTE = {"navy": "#2E5090", "red": "#C0392B", "teal": "#1A8F8F", "gold": "#D4A017"}
plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold"})

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

# Left: Moran scatterplot of residuals
ax = axes[0]
resid = df["residuals"].values
resid_std = (resid - resid.mean()) / resid.std()
lag_resid = w.sparse.dot(resid_std)
ax.scatter(resid_std, lag_resid, alpha=0.6, s=32, color=PALETTE["navy"], edgecolor="white", linewidth=0.3)
z = np.polyfit(resid_std, lag_resid, 1)
xs = np.linspace(resid_std.min(), resid_std.max(), 100)
ax.plot(xs, np.poly1d(z)(xs), color=PALETTE["red"], linewidth=2, linestyle="--")
ax.axhline(0, color="grey", linewidth=0.8)
ax.axvline(0, color="grey", linewidth=0.8)
ax.set_xlabel("OLS residual (standardized)")
ax.set_ylabel("Spatially-lagged residual\n(avg. of 5 nearest neighbors)")
ax.set_title(f"Moran Scatterplot: Residual Clustering\n(Moran's I = {moran_resid.I:.2f}, p = {moran_resid.p_sim:.3f})")

# Right: OLS vs Spatial Lag coefficients
ax = axes[1]
vars_ = ["HDI","HPV_vax_\ncoverage","Screening\nprogram","Smoking\nprev.","HIV\nprev."]
ols_coef = [-4.058, 0.0052, 0.331, 0.0027, 0.0750]
lag_coef = [-2.092, 0.0024, 0.168, 0.0012, 0.0313]
x = np.arange(len(vars_))
width = 0.35
# normalize for visual comparability isn't ideal since different scales; instead show % change
pct_shrink = [100*(1 - abs(l)/abs(o)) for o,l in zip(ols_coef, lag_coef)]
bars = ax.bar(vars_, pct_shrink, color=PALETTE["gold"], edgecolor="white")
ax.set_ylabel("% shrinkage in coefficient magnitude\n(OLS -> Spatial Lag model)")
ax.set_title("How Much Does Accounting for Spatial\nDependence Shrink Each Effect?")
ax.axhline(0, color="black", linewidth=0.8)
for bar, val in zip(bars, pct_shrink):
    ax.annotate(f"{val:.0f}%", (bar.get_x()+bar.get_width()/2, val), xytext=(0,4),
                textcoords="offset points", ha="center", fontsize=9.5, fontweight="bold")

plt.tight_layout()
plt.savefig("fig_spatial.png", dpi=180, bbox_inches="tight")
print("saved")
