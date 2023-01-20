from experiment import RunConfig, experiment

n_list = [10*i for i in range(1, 41)] 
iter_max = 1000
version_list = ["par", "gpu_map", "gpu_mcp"]
frobenius_list = [False, True]

configs = [RunConfig(n, iter_max, 0, num_threads=16, version=version, frobenius=frobenius, warm_up=True) for n in n_list for version in version_list for frobenius in frobenius_list]

df = experiment(configs, "data/ex8.csv")
