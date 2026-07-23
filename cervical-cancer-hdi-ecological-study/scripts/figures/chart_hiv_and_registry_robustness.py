import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df_hiv = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_hiv_sample.csv")
df_reg = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_registry_quality.csv")

PALETTE = {"navy": "#2E5090", "red": "#C0392B", "teal": "#1A8F8F", "gold": "#D4A017"}
plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold"})

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

# Left: HIV vs ASIR scatter, colored by HDI tier
ax = axes[0]
colors = df_hiv["HDI_tier"].map({"Low": PALETTE["red"], "Medium": PALETTE["gold"], "High": PALETTE["teal"], "Very High": PALETTE["navy"]})
ax.scatter(df_hiv["HIV_prev"], df_hiv["ASIR"], c=colors, alpha=0.7, s=35, edgecolor="white", linewidth=0.3)
ax.set_xlabel("HIV prevalence, ages 15-49 (%)")
ax.set_ylabel("ASIR (per 100,000)")
ax.set_title("HIV Prevalence vs. Cervical Cancer Incidence (n=131)")
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=8, label=t)
           for t, c in zip(["Low","Medium","High","Very High"], [PALETTE["red"],PALETTE["gold"],PALETTE["teal"],PALETTE["navy"]])]
ax.legend(handles=handles, fontsize=8.5, loc="upper left", frameon=False, title="HDI tier")

# Right: HDI-ASIR correlation strength by registry quality stratum
ax = axes[1]
strata = ["High (>=90%)", "Lower (<90%)"]
for i, grp in enumerate(strata):
    sub = df_reg[df_reg["registry_quality"]==grp]
    color = PALETTE["navy"] if i==0 else PALETTE["gold"]
    ax.scatter(sub["HDI"], sub["ASIR"], alpha=0.6, s=32, color=color, edgecolor="white", linewidth=0.3, label=f"{grp} (n={len(sub)})")
    z = np.polyfit(sub["HDI"], np.log(sub["ASIR"]), 1)
    xs = np.linspace(sub["HDI"].min(), sub["HDI"].max(), 100)
    ax.plot(xs, np.exp(np.poly1d(z)(xs)), color=color, linewidth=2, linestyle="--")
ax.set_xlabel("Human Development Index")
ax.set_ylabel("ASIR (per 100,000)")
ax.set_title("HDI-Incidence Relationship Holds Regardless\nof Death-Registration Data Quality")
ax.legend(fontsize=8.5, loc="upper right", frameon=False)

plt.tight_layout()
plt.savefig("fig_robustness.png", dpi=180, bbox_inches="tight")
print("saved")
