from experiment import RunConfig, experiment, configs_to_df

num_threads_list = [1] + list(range(2, 17, 2))
parallel_list = [0, 1, 2, 3]
gauss_seidel_list = [True, False]
configs = [RunConfig(gauss_seidel=gauss_seidel, n=200, iter_max=2000, tolerance=0, parallel=parallel, num_threads=num_threads, proc_bind=proc_bind) for proc_bind in ["close", "spread"] for num_threads in num_threads_list for parallel in range(0, 3) for gauss_seidel in [True, False]]


experiment(configs, "data.csv")
