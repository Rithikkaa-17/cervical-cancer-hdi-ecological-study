import pandas as pd
import numpy as np
from scipy import stats

panel = pd.read_csv("cervical_cancer_panel_2022_2024.csv")

def hdi_tier(h):
    if h >= 0.8: return "Very High"
    elif h >= 0.7: return "High"
    elif h >= 0.55: return "Medium"
    else: return "Low"

panel["HDI_tier"] = panel["HDI_val"].apply(hdi_tier)

print("="*65)
print("SPEARMAN CORRELATION: HDI vs ASIR, BY CYCLE")
print("="*65)
for yr in [2022, 2024]:
    sub = panel[panel["Year"]==yr]
    rho, p = stats.spearmanr(sub["HDI_val"], sub["ASIR_val"])
    print(f"{yr} (n={len(sub)}): Spearman rho = {rho:.3f}, p = {p:.2e}")

print("\n" + "="*65)
print("MEDIAN ASIR BY HDI TIER, BY CYCLE")
print("="*65)
tier_order = ["Low","Medium","High","Very High"]
summary = panel.groupby(["Year","HDI_tier"])["ASIR_val"].agg(["median","count"]).reindex(
    pd.MultiIndex.from_product([[2022,2024], tier_order], names=["Year","HDI_tier"]))
print(summary)

print("\n" + "="*65)
print("LOW-HDI vs VERY-HIGH-HDI RATIO, BY CYCLE")
print("="*65)
for yr in [2022, 2024]:
    sub = panel[panel["Year"]==yr]
    low_med = sub[sub["HDI_tier"]=="Low"]["ASIR_val"].median()
    vhigh_med = sub[sub["HDI_tier"]=="Very High"]["ASIR_val"].median()
    print(f"{yr}: Low HDI median={low_med:.1f}, Very High HDI median={vhigh_med:.1f}, Ratio={low_med/vhigh_med:.2f}x")

# --- Countries present in BOTH cycles: paired comparison ---
common = panel.pivot_table(index="Code", columns="Year", values="ASIR_val")
common = common.dropna()
common["change"] = common[2024] - common[2022]
common["pct_change"] = 100 * common["change"] / common[2022]
print(f"\n{'='*65}")
print(f"PAIRED COUNTRIES (n={len(common)}): ASIR CHANGE 2022->2024")
print("="*65)
print(f"Median change: {common['change'].median():+.2f} per 100,000 ({common['pct_change'].median():+.1f}%)")
print(f"Countries with INCREASE: {(common['change']>0).sum()} | DECREASE: {(common['change']<0).sum()}")

# Merge with HDI tier for paired analysis
hdi_lookup = panel[panel["Year"]==2024][["Code","HDI_tier","HDI_val"]].set_index("Code")
common = common.join(hdi_lookup)
print("\nMedian % change in ASIR by HDI tier (paired countries only):")
print(common.groupby("HDI_tier")["pct_change"].median().reindex(tier_order))

common.to_csv("paired_country_changes.csv")
