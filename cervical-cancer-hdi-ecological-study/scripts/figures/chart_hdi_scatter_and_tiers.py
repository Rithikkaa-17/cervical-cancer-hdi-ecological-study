import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_with_HDI.csv")

plt.rcParams["font.family"] = "DejaVu Sans"

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Scatter: HDI vs ASIR
ax = axes[0]
ax.scatter(df["HDI"], df["ASIR"], alpha=0.6, s=35, color="#2E5090", edgecolor="white", linewidth=0.4)
z = np.polyfit(df["HDI"], df["ASIR"], 1)
xs = np.linspace(df["HDI"].min(), df["HDI"].max(), 100)
ax.plot(xs, np.poly1d(z)(xs), color="#C0392B", linewidth=2, linestyle="--")
# label a few extremes
for _, row in df.nlargest(3, "ASIR").iterrows():
    ax.annotate(row["Country"], (row["HDI"], row["ASIR"]), fontsize=8, xytext=(4,2), textcoords="offset points")
ax.set_xlabel("Human Development Index (2023)", fontsize=11)
ax.set_ylabel("Age-Standardized Incidence Rate\n(per 100,000 women-years)", fontsize=11)
ax.set_title("HDI vs. Cervical Cancer Incidence (n=176)", fontsize=12, fontweight="bold")
ax.text(0.05, 0.95, "Spearman ρ = -0.63, p < 0.001", transform=ax.transAxes, fontsize=9.5, va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="#ccc"))
ax.spines[["top","right"]].set_visible(False)

# Bar: median by HDI tier
ax2 = axes[1]
tier_order = ["Low","Medium","High","Very High"]
summary = df.groupby("HDI_tier")[["ASIR","ASMR"]].median().reindex(tier_order)
x = np.arange(len(tier_order))
width = 0.35
b1 = ax2.bar(x - width/2, summary["ASIR"], width, label="Incidence (ASIR)", color="#2E5090")
b2 = ax2.bar(x + width/2, summary["ASMR"], width, label="Mortality (ASMR)", color="#C0392B")
ax2.set_xticks(x)
ax2.set_xticklabels(tier_order)
ax2.set_ylabel("Rate per 100,000 women-years", fontsize=11)
ax2.set_title("Median Burden by HDI Tier", fontsize=12, fontweight="bold")
ax2.legend(frameon=False, fontsize=9.5)
ax2.spines[["top","right"]].set_visible(False)
for bars in [b1,b2]:
    for bar in bars:
        h = bar.get_height()
        ax2.annotate(f"{h:.1f}", (bar.get_x()+bar.get_width()/2, h), textcoords="offset points",
                     xytext=(0,3), ha="center", fontsize=8.5)

plt.tight_layout()
plt.savefig("/home/claude/paper/hdi_cervical_chart.png", dpi=200, bbox_inches="tight")
print("saved")
