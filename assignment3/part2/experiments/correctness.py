from experiment import RunConfig, experiment

n_list = [30, 100]
iter_max = 10000
version_list = ["seq", "par"]

configs = [RunConfig(n, iter_max, 0, test=True, num_threads=16) for n in n_list for version in version_list]

df = experiment(configs, "data/data.csv")
print(f"Maximum error: {df['error'].max()}")
