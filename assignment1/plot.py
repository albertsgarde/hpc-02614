import sys
import os
import math
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path

if len(sys.argv) < 2:
    raise Exception("Script takes two arguments: plot.py <data_path> <out_path>")

data_path = sys.argv[1]
perf_data_path = os.path.join(data_path, "perf.dat")
cache_data_path = os.path.join(data_path, "cache.dat")
all_data_path = os.path.join(data_path, "data.csv")
out_path = sys.argv[2]

Path(out_path).mkdir(parents=True, exist_ok=True)


perf_df = pd.read_csv(perf_data_path, delimiter=" ", names=["Memory", "MFLOPs", "Checksum", "Hash", "Type", "SizeOrBlock", "SizeOrNA"], skipinitialspace=True, index_col=False)
perf_df[["Block", "Size"]] = perf_df[["SizeOrBlock", "SizeOrNA"]].apply(lambda row: [1, row[0]] if math.isnan(row[1]) else row, axis=1, raw=True)
perf_df["Type"] = perf_df["Type"].apply(lambda type: type[8:])
perf_df = perf_df.drop(columns=["Hash", "SizeOrBlock", "SizeOrNA"])



cache_df = pd.read_csv(cache_data_path, header=0, names=["Time exclusive", "Time", "L1H", "L1M", "L2H", "L2M", "L3H", "L3M", "Function", "Type", "Size", "Block"], index_col=False)
cache_df = cache_df.drop(columns=["Time exclusive", "Function"])

df = perf_df.merge(cache_df, on=["Type", "Size", "Block"])

for cache in range(1, 4):
    df[f"L{cache}R"] = df[[f"L{cache}H", f"L{cache}M"]].apply(lambda row: row[1]/row[0] if row[0] != 0 else 0, axis=1)

df = df[["Type", "Size", "Block", "Checksum", "Memory", "Time", "MFLOPs", "L1H", "L1M", "L1R", "L2H", "L2M", "L2R", "L3H", "L3M", "L3R"]]

df.to_csv(all_data_path)

x_axes = ["Size", "Memory"]
measurements = ["Time", "MFLOPs", "L1R", "L2R", "L3R"]

for x_axis in x_axes:
    for measurement in measurements:
        plot_path = os.path.join(out_path, f"{measurement}_{x_axis}.pdf")
        #["Type", "Block"]
        plot = sns.lineplot(df, x=x_axis, y=measurement, hue="Type")
        plt.xscale('log')
        plot.get_figure().savefig(plot_path)
