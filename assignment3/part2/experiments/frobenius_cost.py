from experiment import RunConfig, experiment

n_list = [100, 200] 
iter_max = 10000
version_list = ["par", "gpu_mcp"]
frobenius_list = [True, False]

configs = [RunConfig(n, iter_max, 0, num_threads=16, version=version, frobenius=frobenius) for n in n_list for version in version_list for frobenius in frobenius_list]

df = experiment(configs, "data/data.csv")
