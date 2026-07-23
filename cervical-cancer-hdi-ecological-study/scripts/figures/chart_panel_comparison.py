import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

panel = pd.read_csv("cervical_cancer_panel_2022_2024.csv")
common = pd.read_csv("paired_country_changes.csv")

PALETTE = {"navy": "#2E5090", "red": "#C0392B", "teal": "#1A8F8F", "gold": "#D4A017", "grey": "#7f8c8d"}
plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold"})

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

# Left: paired scatter, 2022 vs 2024
ax = axes[0]
colors = common["HDI_tier"].map({"Low": PALETTE["red"], "Medium": PALETTE["gold"], "High": PALETTE["teal"], "Very High": PALETTE["navy"]})
ax.scatter(common["2022"], common["2024"], c=colors, alpha=0.7, s=35, edgecolor="white", linewidth=0.3)
lims = [0, max(common["2022"].max(), common["2024"].max())+5]
ax.plot(lims, lims, color="black", linestyle="--", linewidth=1.2, label="No change (y=x)")
ax.set_xlim(lims); ax.set_ylim(lims)
ax.set_xlabel("ASIR 2022 (per 100,000)")
ax.set_ylabel("ASIR 2024 (per 100,000)")
ax.set_title("Same-Country Incidence: 2022 vs. 2024 (n=156)")
handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=8, label=t)
           for t, c in zip(["Low","Medium","High","Very High"], [PALETTE["red"],PALETTE["gold"],PALETTE["teal"],PALETTE["navy"]])]
ax.legend(handles=handles+[plt.Line2D([0],[0], color='black', linestyle='--', label='No change')], fontsize=8.5, loc="upper left", frameon=False)

# Right: % change by HDI tier
ax = axes[1]
tier_order = ["Low","Medium","High","Very High"]
data = [common.loc[common["HDI_tier"]==t, "pct_change"].dropna() for t in tier_order]
bp = ax.boxplot(data, labels=tier_order, patch_artist=True, widths=0.55, showfliers=False)
for patch, c in zip(bp["boxes"], [PALETTE["red"],PALETTE["gold"],PALETTE["teal"],PALETTE["navy"]]):
    patch.set_facecolor(c); patch.set_alpha(0.65)
ax.axhline(0, color="black", linewidth=1, linestyle="-")
ax.set_ylabel("% change in ASIR, 2022 -> 2024")
ax.set_title("Incidence Change by HDI Tier (Paired Countries)")
plt.setp(ax.get_xticklabels(), rotation=10)

plt.tight_layout()
plt.savefig("fig_panel_2022_2024.png", dpi=180, bbox_inches="tight")
print("saved")

# Print summary stats for text
for t in tier_order:
    vals = common.loc[common["HDI_tier"]==t, "pct_change"].dropna()
    print(f"{t}: median={vals.median():+.2f}%, n={len(vals)}")
