import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors

world = gpd.read_file("countries.geojson")
print(world.shape, world.columns.tolist())

df = pd.read_csv("/mnt/user-data/outputs/cervical_cancer_full_covariates.csv")
merged = world.merge(df[["code","ASIR"]], left_on="ISO3166-1-Alpha-3", right_on="code", how="left")
print("Matched:", merged["ASIR"].notna().sum(), "/", len(df))

fig, ax = plt.subplots(figsize=(15, 8))
cmap = LinearSegmentedColormap.from_list("cervix", ["#e8f0fb","#a8c5e8","#5b8fc9","#2E5090","#c0392b","#7b1c14"])

merged.plot(column="ASIR", cmap=cmap, linewidth=0.3, ax=ax, edgecolor="#888888",
            missing_kwds={"color": "#e8e8e8", "edgecolor": "#bbbbbb", "linewidth": 0.3, "label": "No data"},
            vmin=0, vmax=95, legend=True,
            legend_kwds={"label": "ASIR (per 100,000 women-years)", "orientation": "horizontal",
                         "shrink": 0.4, "pad": 0.02, "aspect": 30})

ax.set_title("Global Cervical Cancer Incidence, 2024\n(Age-Standardized Incidence Rate per 100,000 Women-Years)",
              fontsize=15, fontweight="bold", pad=15)
ax.set_axis_off()
ax.set_xlim(-170, 190)
ax.set_ylim(-58, 85)

plt.tight_layout()
plt.savefig("fig_world_map.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved")
