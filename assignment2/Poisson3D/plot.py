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
is_gauss_seidel = df['gauss_seidel']
df.loc[df['gauss_seidel'] == True, 'gauss_seidel'] = "gauss_seidel"
df.loc[df['gauss_seidel'] == False, 'gauss_seidel'] = "jacobi"

base_time_gs = df.loc[(df['gauss_seidel'] == "gauss_seidel") & (df['num_threads'] == 1), "time"].values[0]
base_time_j = df.loc[(df['gauss_seidel'] == "jacobi") & (df['num_threads'] == 1), "time"].values[0]

df.loc[df['gauss_seidel'] == "gauss_seidel", "speedup"] = base_time_gs/df.loc[:,"time"]
df.loc[df['gauss_seidel'] == "jacobi", "speedup"] = base_time_j/df.loc[:,"time"]

plot_path = os.path.join(out_path, "out.png")

fig = plt.figure()
ax = fig.add_subplot(111)
plot = sns.lineplot(data=df, x="num_threads", y="speedup", hue="gauss_seidel", style="gauss_seidel", markers=True, dashes=False, ax=ax)

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.legend(title="")
plot.get_figure().savefig(plot_path)

