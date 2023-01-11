from experiment import RunConfig, experiment

n_list = range(1, 31)
iter_max = 1000
tolerance = 1e-5

configs = [RunConfig(False, n, iter_max, tolerance) for n in n_list]
configs += [RunConfig(True, n, iter_max, tolerance) for n in n_list]

experiment(configs, "data.csv")
