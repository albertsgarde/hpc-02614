from experiment import RunConfig, experiment

n_list = range(1, 31)
iter_max = 1000
tolerance = 1e-5

configs = [RunConfig(False, n, iter_max, tolerance, test=True) for n in n_list]
configs += [RunConfig(True, n, iter_max, tolerance, test=True) for n in n_list]

df = experiment(configs, "data.csv")
print(f"Maximum error: {df['error'].max()}")
