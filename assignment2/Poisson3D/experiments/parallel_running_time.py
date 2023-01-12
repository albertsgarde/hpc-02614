from experiment import RunConfig, experiment

num_threads_list = [1] + list(range(2, 17, 2))
configs = [RunConfig(gauss_seidel=gauss_seidel, n=100, iter_max=10000, tolerance=0, parallel=parallel, num_threads=num_threads) for num_threads in num_threads_list for parallel in range(0, 3) for gauss_seidel in [True, False]]


experiment(configs, "data.csv")
