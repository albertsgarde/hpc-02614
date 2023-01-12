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
import numpy as np

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

f = 0.99
num_threads = df.num_threads.unique()
amdahl = 1.0/((1.0 - f) + f/num_threads)
sns.lineplot(x=num_threads, y=amdahl, color="orange", label=f"amdahl f={f}")
f = 0.9
num_threads = df.num_threads.unique()
amdahl = 1.0/((1.0 - f) + f/num_threads)
sns.lineplot(x=num_threads, y=amdahl, color="blue", label=f"amdahl f={f}")


sns.lineplot(data=df, x="num_threads", y="num_threads", color="green", label="linear",linestyle='--', ax=ax)
plot = sns.lineplot(data=df, x="num_threads", y="speedup", hue="gauss_seidel", style="gauss_seidel", linestyle='', marker='o', dashes=False, ax=ax)

ax.yaxis.set_major_formatter(ScalarFormatter())
ax.legend(title="")
ax.set_xlabel("Num. threads")
ax.set_ylabel("Speedup")
ax.set_title("Speedup for OMP_PROC_BIND=spread and N=100")
plot.get_figure().savefig(plot_path)