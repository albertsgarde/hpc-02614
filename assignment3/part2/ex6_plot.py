import sys
import os
import math
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import SymLogNorm
import pandas as pd
import seaborn as sns
from pathlib import Path

df = pd.read_csv("data/ex6.csv")

def host_row(df, row):
    par1 = df[(df["version"] == "par") & (df["n"] == row["n"]) & (df["frobenius"] == row["frobenius"]) & (df["iter_max"] == row["iter_max"]) & (df["warm_up"] == row["warm_up"])]
    if len(par1) > 1:
        raise Exception("Multiple rows found")
    if len(par1) < 1:
        raise Exception("No rows found")
    return par1.iloc[0]

df["speedup"] = df.apply(lambda row: row["iterations_per_second"]/host_row(df, row)["iterations_per_second"], axis=1)

plot_path = os.path.join("plots", "ex6.pdf")


fig = plt.figure()
ax = fig.add_subplot(111)
plt.title("Speedup")
plot = sns.lineplot(data=df, x="n", y="speedup", hue="version", markers=True, dashes=False, ax=ax)
ax.yaxis.set_major_formatter(ScalarFormatter())
plot.get_figure().savefig(plot_path)