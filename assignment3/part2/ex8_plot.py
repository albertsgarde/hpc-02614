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

df = pd.read_csv("data/ex8.csv")

# Speedup
def host_row(df, row):
    par1 = df[(df["version"] == "par") & (df["n"] == row["n"]) & (df["frobenius"] == row["frobenius"]) & (df["iter_max"] == row["iter_max"]) & (df["warm_up"] == row["warm_up"])]
    if len(par1) > 1:
        raise Exception("Multiple rows found")
    if len(par1) < 1:
        raise Exception("No rows found")
    return par1.iloc[0]

df["speedup"] = df.apply(lambda row: row["iterations_per_second"]/host_row(df, row)["iterations_per_second"], axis=1)

plot_path = os.path.join("plots", "ex8_speedup.pdf")


fig = plt.figure()
ax = fig.add_subplot(111)
plt.title("Speedup")
plot = sns.lineplot(data=df[df["frobenius"] == True], x="n", y="speedup", hue="version", markers=True, dashes=False, ax=ax)
ax.yaxis.set_major_formatter(ScalarFormatter())
plot.get_figure().savefig(plot_path)


# Frobenius cost
cost_df = df[['n', 'frobenius', 'version', 'iterations_per_second']]
cost_pivot_df = cost_df.pivot(index=["n", "version"], columns="frobenius", values="iterations_per_second")


test_df = cost_df
test_pivot_df = test_df.pivot(index=["n", "version"], columns="frobenius", values="iterations_per_second")
slowdown_series = test_pivot_df[False]/test_pivot_df[True]
slowdown_df = slowdown_series.to_frame().reset_index().rename(columns= {0: 'slowdown'})

plot_path = os.path.join("plots", "ex8_frobenius_cost.pdf")

fig = plt.figure()
ax = fig.add_subplot(111)
plt.title("Cost of Frobenius norm")
plot = sns.lineplot(data=slowdown_df, x="n", y="slowdown", hue="version", markers=True, dashes=False, ax=ax)
ax.yaxis.set_major_formatter(ScalarFormatter())
plot.get_figure().savefig(plot_path)