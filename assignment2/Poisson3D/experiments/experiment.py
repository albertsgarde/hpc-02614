import sys
import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
import subprocess

class RunConfig:
    def __init__(self, gauss_seidel: bool, n: int, iter_max: int, tolerance: float = 0, test: bool = False, parallel: int = 0, num_threads: int = 1, schedule: str = "", proc_bind: str = ""):
        self.gauss_seidel = gauss_seidel
        self.n = n
        self.iter_max = iter_max
        self.tolerance = tolerance
        self.test = test
        self.parallel = parallel
        self.num_threads = num_threads
        self.schedule = schedule
        self.proc_bind = proc_bind
    
    def to_row(self):
        return [self.gauss_seidel, self.n, self.iter_max, self.tolerance, self.test, self.parallel, self.num_threads, self.schedule, self.proc_bind]

def run(config: RunConfig):
    env_vars = { "OMP_NUM_THREADS": str(config.num_threads) }
    if config.schedule != "":
        env_vars["OMP_SCHEDULE"] = config.schedule
    if config.proc_bind != "":
        env_vars["OMP_PROC_BIND"] = config.proc_bind
    env = {**os.environ.copy(), **env_vars}
    bin = "poisson_" + ("gs" if config.gauss_seidel else "j")
    args = [bin, f"{config.n}", f"{config.iter_max}", f"{config.tolerance}", "0"] + (["-t"] if config.test else []) + ["-p", f"{config.parallel}"]
    print(f"Running `{' '.join([f'{key}={value}' for key,value in env_vars.items()])} {' '.join(args)}`")
    process = subprocess.Popen(args, env=env, stdout=subprocess.PIPE)

    out, err = process.communicate()
    if process.returncode != 0:
        print(err)
        sys.exit(1)
    print(str(out))
    [time_string, num_iterations_string, error_string] = out.split(b' ')
    time = float(time_string)
    num_iterations = int(num_iterations_string)
    error = float(error_string)
    row = config.to_row() + [time, num_iterations, error]
    return row

def experiment(configs, data_path: str):
    print(f"Running {len(configs)} experiments...")
    rows = [run(config) for config in configs]
    df = pd.DataFrame(rows, columns=["gauss_seidel", "n", "iter_max", "tolerance", "test", "parallel", "num_threads", "schedule", "proc_bind", "time", "num_iterations", "error"])
    df["iterations_per_second"] = df[["time", "num_iterations"]].apply(lambda row: row[1]/row[0], axis=1)
    df.to_csv(data_path, index=False)
    return df

def configs_to_df(configs):
    rows = [config.to_row() for config in configs]
    df = pd.DataFrame(rows, columns=["gauss_seidel", "n", "iter_max", "tolerance", "test", "parallel", "num_threads", "schedule"])
    return df 
