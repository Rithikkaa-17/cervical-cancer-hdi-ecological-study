"""
Reproduce every statistic reported in the manuscript
"Socioeconomic, Health-System, and Spatial Determinants of Cervical Cancer Burden:
A Global Ecological Study" from the processed datasets in data/processed/.

Usage (from the repository root):
    pip install pandas numpy scipy statsmodels libpysal esda spreg geopandas
    python scripts/reproduce_all.py

Queen-contiguity weights need country boundary polygons; the script downloads
data/raw/countries.geojson (datasets/geo-countries) if it is not already present.
"""
import os, json, warnings, urllib.request
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy import stats
warnings.filterwarnings("ignore")

P = "data/processed/"
out = {}
def show(title): print("\n" + "=" * 70 + "\n" + title + "\n" + "=" * 70)

# ---------------------------------------------------------------- Section IV
full = pd.read_csv(P + "02_with_hdi.csv")                    # n = 176
show(f"IV-A  Descriptives (n = {len(full)})")
for v in ["ASIR", "ASMR", "MIR"]:
    q = full[v].quantile([.25, .5, .75])
    print(f"{v}: median {q[.5]:.2f}  IQR {q[.25]:.2f}-{q[.75]:.2f}  range {full[v].min()}-{full[v].max()}")

show("IV-B  Table III: HDI vs outcomes (n = 176)")
for v in ["ASIR", "ASMR", "MIR"]:
    rho, p = stats.spearmanr(full.HDI, full[v]); r, pr = stats.pearsonr(full.HDI, full[v])
    line = f"{v}: Spearman {rho:.2f} (p={p:.1e})  Pearson {r:.2f} (p={pr:.1e})"
    if v != "MIR":
        m = smf.ols(f"np.log({v}) ~ HDI", full).fit()
        line += f"  log-OLS HDI {m.params.HDI:.2f} (p={m.pvalues.HDI:.1e}) R2 {m.rsquared:.3f}"
    print(line)
med = full.groupby("HDI_tier")[["ASIR", "ASMR"]].median()
print(med.round(2))
print(f"Low:VeryHigh ratio  ASIR {med.ASIR['Low']/med.ASIR['Very High']:.2f}x  ASMR {med.ASMR['Low']/med.ASMR['Very High']:.2f}x")

cov = pd.read_csv(P + "03_full_covariates_n156.csv")         # n = 156
X = "HDI + HPV_vax_coverage + Screening_program_bin + Smoking_prev_female"
show(f"IV-D  Table IV: multivariable OLS (n = {len(cov)})")
for y in ["logASIR", "logASMR"]:
    m = smf.ols(f"{y} ~ {X}", cov).fit()
    print(y, " ".join(f"{k}={m.params[k]:+.4f}(p={m.pvalues[k]:.3f})" for k in m.params.index[1:]), f"R2={m.rsquared:.2f}")
Xm = cov[["HDI", "HPV_vax_coverage", "Screening_program_bin", "Smoking_prev_female"]].assign(c=1)
print("VIF:", {c: round(variance_inflation_factor(Xm.values, i), 2) for i, c in enumerate(Xm.columns[:-1])})
for v in ["logASIR", "logASMR"]:
    rho, p = stats.spearmanr(cov.Smoking_prev_female, cov[v]); print(f"smoking vs {v}: rho {rho:.2f} p {p:.1e}")

# ---------------------------------------------------------------- Section V-A
show("V-A  Temporal comparison")
panel = pd.read_csv(P + "06_panel_2022_2024.csv")
def tier(h): return pd.cut(h, [0, .55, .7, .8, 1.01], right=False, labels=["Low", "Medium", "High", "Very High"])
for yr in (2022, 2024):
    s = panel[panel.Year == yr]; m = s.groupby(tier(s.HDI_val), observed=True).ASIR_val.median()
    print(f"{yr} (n={len(s)}): naive Low:VeryHigh ASIR ratio {m['Low']/m['Very High']:.3f}x")
pc = pd.read_csv(P + "07_paired_country_changes.csv")
print(f"paired n={len(pc)} median %chg {pc['pct_change'].median():.2f} up {(pc['pct_change']>0).sum()} down {(pc['pct_change']<0).sum()}")
print("Wilcoxon p", round(stats.wilcoxon(pc["2024"], pc["2022"]).pvalue, 3),
      " Kruskal-Wallis", stats.kruskal(*[g.values for _, g in pc.groupby("HDI_tier")["pct_change"]]))
print("tier medians %chg", pc.groupby("HDI_tier")["pct_change"].median().round(2).to_dict())
rho, p = stats.spearmanr(pc.HDI_val, pc['pct_change']); print(f"HDI vs %chg rho {rho:.3f} p {p:.3f}")
g22 = panel[panel.Year == 2022].set_index("Code")
both = full.set_index("Alpha-3 code").join(g22[["ASIR_val"]], how="inner")
chg = (both.ASIR - both.ASIR_val) / both.ASIR_val * 100
print(f"full-sample paired n={len(both)} median {chg.median():.2f}% Wilcoxon p {stats.wilcoxon(both.ASIR, both.ASIR_val).pvalue:.3f}"
      f" KW p {stats.kruskal(*[g.values for _, g in chg.groupby(both.HDI_tier)]).pvalue:.3f}")

# ---------------------------------------------------------------- Section V-B
hiv = pd.read_csv(P + "04_hiv_mediation_sample_n131.csv")    # n = 131
show(f"V-B  HIV-extended model (n = {len(hiv)})")
for v in ["ASIR", "ASMR"]:
    rho, p = stats.spearmanr(hiv.HIV_prev, hiv[v]); print(f"HIV vs {v}: rho {rho:.2f} p {p:.1e}")
for y in ["logASIR", "logASMR"]:
    m0 = smf.ols(f"{y} ~ {X}", hiv).fit(); m1 = smf.ols(f"{y} ~ {X} + HIV_prev", hiv).fit()
    print(y, f"HIV b={m1.params.HIV_prev:.3f} p={m1.pvalues.HIV_prev:.1e} R2 {m0.rsquared:.2f}->{m1.rsquared:.2f} "
          f"HDI {m0.params.HDI:.2f}->{m1.params.HDI:.2f} (attenuation {100*(1-m1.params.HDI/m0.params.HDI):.1f}%)",
          " others:", {k: round(m1.pvalues[k], 3) for k in ["HPV_vax_coverage", "Screening_program_bin"]})
Xh = hiv[["HDI", "HPV_vax_coverage", "Screening_program_bin", "Smoking_prev_female", "HIV_prev"]].assign(c=1)
print("VIF HIV:", round(variance_inflation_factor(Xh.values, 4), 2))

reg = pd.read_csv(P + "05_registry_quality_sample_n142.csv")
print(f"\nRegistry strata (n={len(reg)}):", reg.registry_quality.value_counts().to_dict())
rho, p = stats.spearmanr(reg.HDI, reg.Death_reg_pct); print(f"HDI vs registration rho {rho:.2f} p {p:.1e}")
for g, s in reg.groupby("registry_quality"):
    print(g, "HDI-ASIR", round(stats.spearmanr(s.HDI, s.ASIR)[0], 2), "HDI-ASMR", round(stats.spearmanr(s.HDI, s.ASMR)[0], 2),
          "tiers", s.HDI_tier.value_counts().to_dict())

# ---------------------------------------------------------------- Section V-C
show("V-C  Mediation (consistent covariate adjustment)")
Xc = "HPV_vax_coverage + Screening_program_bin + Smoking_prev_female"
ma = smf.ols(f"HIV_prev ~ HDI + {Xc}", hiv).fit(); mb = smf.ols(f"logASMR ~ HDI + HIV_prev + {Xc}", hiv).fit()
mc = smf.ols(f"logASMR ~ HDI + {Xc}", hiv).fit()
a, b = ma.params.HDI, mb.params.HIV_prev; ab = a * b
se = np.sqrt(b**2 * ma.bse.HDI**2 + a**2 * mb.bse.HIV_prev**2)
print(f"a={a:.2f} (p={ma.pvalues.HDI:.3f}) b={b:.4f} ab={ab:.3f} Sobel z={ab/se:.2f} p={2*stats.norm.sf(abs(ab/se)):.3f}")
print(f"c={mc.params.HDI:.2f} c'={mb.params.HDI:.2f} proportion={100*ab/mc.params.HDI:.1f}%")
C = hiv[["HPV_vax_coverage", "Screening_program_bin", "Smoking_prev_female"]].values
X1 = np.column_stack([np.ones(len(hiv)), hiv.HDI, C]); X2 = np.column_stack([X1, hiv.HIV_prev])
def ab_i(i): return (np.linalg.lstsq(X1[i], hiv.HIV_prev.values[i], rcond=None)[0][1] *
                     np.linalg.lstsq(X2[i], hiv.logASMR.values[i], rcond=None)[0][-1])
rng = np.random.default_rng(2026); n = len(hiv)
bs = np.array([ab_i(rng.integers(0, n, n)) for _ in range(5000)])
z0 = stats.norm.ppf((bs < ab).mean()); lo, hi = stats.norm.cdf(2 * z0 + stats.norm.ppf([.025, .975]))
print("bias-corrected bootstrap 95% CI", np.quantile(bs, [lo, hi]).round(3))

# ---------------------------------------------------------------- Section V-D
show("V-D  Spatial analysis (n = 131)")
import libpysal, esda
from spreg import OLS, ML_Lag
sp = pd.read_csv(P + "08_spatial_analysis_sample.csv").reset_index(drop=True)
coords = np.column_stack([sp.longitude, sp.latitude])
yv = sp[["logASIR"]].values
xcols = ["HDI", "HPV_vax_coverage", "Screening_program_bin", "Smoking_prev_female", "HIV_prev"]
xv = sp[xcols].values

def run(w, label):
    w.transform = "r"
    ols = OLS(yv, xv, w=w, spat_diag=True, name_x=xcols, name_y="logASIR")
    mi = esda.Moran(ols.u.flatten(), w, permutations=999)
    lag = ML_Lag(yv, xv, w=w, name_x=xcols, name_y="logASIR")
    names = ["CONST"] + xcols + ["W_logASIR"]
    zs = dict(zip(names, lag.z_stat))
    hdi = lag.betas[1][0]; rho = lag.betas[-1][0]
    print(f"{label:28s} Moran(resid) {mi.I:.3f} p={mi.p_sim:.3f} | lag HDI {hdi:.2f} p={zs['HDI'][1]:.1e} "
          f"rho {rho:.3f} p={zs['W_logASIR'][1]:.1e} pseudoR2 {lag.pr2:.2f} | OLS R2 {ols.r2:.2f} HDI {ols.betas[1][0]:.2f}"
          f" | vax p={zs['HPV_vax_coverage'][1]:.3f} screen p={zs['Screening_program_bin'][1]:.3f} HIV p={zs['HIV_prev'][1]:.3f}"
          f" | total HDI {hdi/(1-rho):.2f}")
    return ols

w5 = libpysal.weights.KNN.from_array(coords, k=5)
ols = run(w5, "k-NN k=5 (primary)")
w5.transform = "r"
print("Moran raw ASIR", round(esda.Moran(sp.ASIR.values, w5, permutations=999).I, 2),
      " Moran HDI", round(esda.Moran(sp.HDI.values, w5, permutations=999).I, 2))
print("robust LM-lag", np.round(ols.rlm_lag, 3), " robust LM-error", np.round(ols.rlm_error, 3))
for k in (3, 8, 10): run(libpysal.weights.KNN.from_array(coords, k=k), f"k-NN k={k}")

# great-circle distances for distance-based weights
R = 6371.0
lat, lon = np.radians(sp.latitude.values), np.radians(sp.longitude.values)
dlat = lat[:, None] - lat[None, :]; dlon = lon[:, None] - lon[None, :]
D = 2 * R * np.arcsin(np.sqrt(np.sin(dlat/2)**2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon/2)**2))
np.fill_diagonal(D, np.inf)
cut = D.min(1).max(); print(f"\nsmallest cutoff giving every country a neighbor: {cut:.0f} km")
def dist_w(c, inv):
    nb, wt = {}, {}
    for i in range(len(D)):
        j = np.where(D[i] <= c)[0]; nb[i] = list(j); wt[i] = list(1 / D[i, j]) if inv else [1.0] * len(j)
    return libpysal.weights.W(nb, wt, silence_warnings=True)
for c in (cut, 3700):
    run(dist_w(c, True), f"inverse-distance {c:.0f} km"); run(dist_w(c, False), f"distance-band {c:.0f} km")

gj = "data/raw/countries.geojson"
try:
    import geopandas as gpd
    if not os.path.exists(gj):
        urllib.request.urlretrieve("https://raw.githubusercontent.com/datasets/geo-countries/master/data/countries.geojson", gj)
    g = gpd.read_file(gj)
    iso = [c for c in g.columns if c.upper() in ("ISO3166-1-ALPHA-3", "ISO_A3", "ADM0_A3")][0]
    g = g[[iso, "geometry"]].rename(columns={iso: "code"}).dissolve("code").reset_index()
    sub = sp.merge(g, on="code", how="inner"); print(f"\ncountries with boundary polygons: {len(sub)}")
    sub = gpd.GeoDataFrame(sub, geometry="geometry")
    wq = libpysal.weights.Queen.from_dataframe(sub, use_index=False, silence_warnings=True)
    print("islands (no land neighbor):", len(wq.islands), f"({100*len(wq.islands)/len(sub):.1f}%)")
    keep = [i for i in range(len(sub)) if i not in set(wq.islands)]
    yv_b, xv_b = yv, xv
    yv, xv = sub[["logASIR"]].values[keep], sub[xcols].values[keep]
    subk = sub.iloc[keep].reset_index(drop=True)
    wq2 = libpysal.weights.Queen.from_dataframe(subk, use_index=False, silence_warnings=True)
    run(wq2, f"Queen contiguity (n={len(subk)}, islands dropped)")
    yv, xv = yv_b, xv_b
except Exception as e:
    print("Queen contiguity skipped:", e)
