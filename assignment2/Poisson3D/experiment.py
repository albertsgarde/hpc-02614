import sys
import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import subprocess

if len(sys.argv) < 2:
    raise Exception("Script takes one argument: experiment.py <out_path>")

data_path = sys.argv[1]

def run(test: bool, gauss_seidel: bool, n: int, iter_max: int, tolerance: float):
    env = os.environ.copy()
    bin = "poisson_" + ("gs" if gauss_seidel else "j")
    args = [bin, f"{n}", f"{iter_max}", f"{tolerance}", "0"] + (["-t"] if test else [])
    print(f"Running `{' '.join(args)}`")
    process = subprocess.Popen(args, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    out, err = process.communicate()
    if process.returncode != 0:
        print(err)
        sys.exit(1)
    [time_string, num_iterations_string] = out.split(b' ')
    time = float(time_string)
    num_iterations = int(num_iterations_string)
    print(err)
    row = [test, gauss_seidel, n, iter_max, tolerance, time, num_iterations]
    return row

n_list = range(1, 31)
iter_max = 1000
tolerance = 1e-5

rows = [run(False, False, n, iter_max, tolerance) for n in n_list]
rows += [run(False, True, n, iter_max, tolerance) for n in n_list]

df = pd.DataFrame(rows, columns=["test", "gauss_seidel", "n", "iter_max", "tolerance", "time", "num_iterations"])
df["iterations_per_second"] = df[["time", "num_iterations"]].apply(lambda row: row[1]/row[0], axis=1)

df.to_csv(data_path)
