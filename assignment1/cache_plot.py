import sys
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

if len(sys.argv) < 1:
    raise Exception("Output path must be specified as the first argument.")

df = pd.read_csv("cache.dat", header=0, names=["Time exclusive", "Time", "L1 Cache Misses", "Instructions", "L2 Cache Misses", "L3 Cache Misses", "Type", "Size"])
df["Type"] = df["Type"].apply(lambda type: type[8:])
df.drop(columns=["Time exclusive"])

measurements = ["Time", "L1 Cache Misses", "L2 Cache Misses", "L3 Cache Misses", "Instructions"]

fig, axs = plt.subplots(nrows=5, figsize=(6, 16))
for index, measurement in enumerate(measurements):
    sns.lineplot(df, x="Size", y=measurement, hue="Type", ax=axs[index])
    
plt.savefig(sys.argv[1])
