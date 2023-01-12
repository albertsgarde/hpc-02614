from experiment import RunConfig, experiment

num_threads=[1,2,4,6,8,10,12,16]

configs = [RunConfig(gauss_seidel=True, n=200, iter_max=1000, tolerance=0, num_threads=n) for n in num_threads]
configs += [RunConfig(gauss_seidel=False, n=200, iter_max=1000, tolerance=0, num_threads=n) for n in num_threads]

experiment(configs, "data.csv")