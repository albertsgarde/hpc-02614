from experiment import RunConfig, experiment

configs = [RunConfig(gauss_seidel=True, n=100, iter_max=iter_max, tolerance=0) for iter_max in range(100, 2000, 100)]

experiment(configs, "data.csv")
