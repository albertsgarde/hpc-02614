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


plot_path = os.path.join(out_path, "out.pdf")

fig = plt.figure()
ax = fig.add_subplot(111)
plot = sns.lineplot(df, x="n", y="num_iterations", hue="gauss_seidel", style="gauss_seidel", markers=True, dashes=False, ax=ax)
ax.yaxis.set_major_formatter(ScalarFormatter())
plot.get_figure().savefig(plot_path)
