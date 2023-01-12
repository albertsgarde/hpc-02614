from experiment import RunConfig, experiment

num_threads=[1,2,4,8,16]

configs = [RunConfig(gauss_seidel=True, n=100, iter_max=10000, tolerance=0, num_threads=n) for n in num_threads]
configs += [RunConfig(gauss_seidel=False, n=100, iter_max=10000, tolerance=0, num_threads=n) for n in num_threads]

experiment(configs, "data.csv")
