import sys
import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import subprocess

class RunConfig:
    def __init__(self, gauss_seidel: bool, n: int, iter_max: int, tolerance: float = 0, test: bool = False, parallel: bool = False, num_threads: int = 1, schedule: str = "static"):
        self.gauss_seidel = gauss_seidel
        self.n = n
        self.iter_max = iter_max
        self.tolerance = tolerance
        self.test = test
        self.parallel = parallel
        self.num_threads = num_threads
        self.schedule = schedule

def run(config: RunConfig):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(config.num_threads)
    env["OMP_SCHEDULE"] = config.schedule
    bin = "poisson_" + ("gs" if config.gauss_seidel else "j")
    args = [bin, f"{config.n}", f"{config.iter_max}", f"{config.tolerance}", "0"] + (["-t"] if config.test else [])
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
    row = [config.gauss_seidel, config.n, config.iter_max, config.tolerance, config.test, config.parallel, config.num_threads, config.schedule, time, num_iterations]
    return row

def experiment(configs: list[RunConfig], data_path: str):
    rows = [run(config) for config in configs]
    df = pd.DataFrame(rows, columns=["gauss_seidel", "n", "iter_max", "tolerance", "test", "parallel", "num_threads", "schedule", "time", "num_iterations"])
    df["iterations_per_second"] = df[["time", "num_iterations"]].apply(lambda row: row[1]/row[0], axis=1)
    df.to_csv(data_path)
