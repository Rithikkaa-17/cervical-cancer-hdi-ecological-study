# Cervical Cancer HDI Ecological Study — project files

This folder contains all data, scripts, figures and the notebook for the study.

The full documentation (data sources, pipeline, reproduction instructions and key findings) is in the **[repository README](../README.md)**.

Quick reproduction of every statistic in the manuscript:

```bash
pip install pandas numpy scipy statsmodels libpysal esda spreg geopandas
python scripts/reproduce_all.py
```

Expected output: [`scripts/reproduce_all_output.txt`](scripts/reproduce_all_output.txt).
