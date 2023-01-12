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

if len(sys.argv) < 3:
    raise Exception("Script takes two arguments: plot.py <data_path> <out_path>")

data_path = sys.argv[1]
out_path = sys.argv[2]

df = pd.read_csv(data_path)

def ser_row(df, row):
    par1 = df[(df["gauss_seidel"] == row["gauss_seidel"]) & (df["parallel"] == 0) & (df["num_threads"] == row["num_threads"])]
    if len(par1) > 1:
        raise Exception("Multiple rows found")
    return par1.iloc[0]

df = df[["num_threads", "gauss_seidel", "parallel", "iterations_per_second"]]
df.loc[df['gauss_seidel'] == True, 'gauss_seidel'] = "gauss_seidel"
df.loc[df['gauss_seidel'] == False, 'gauss_seidel'] = "jacobi"

df["speedup"] = df.apply(lambda row: row["iterations_per_second"]/ser_row(df, row)["iterations_per_second"], axis=1)

plot_path = os.path.join(out_path, "out.pdf")



fig = plt.figure()
ax = fig.add_subplot(111)
plot = sns.lineplot(data=df, x="num_threads", y="speedup", hue="parallel", style="gauss_seidel", markers=True, dashes=False, ax=ax)
ax.yaxis.set_major_formatter(ScalarFormatter())
plot.get_figure().savefig(plot_path)
