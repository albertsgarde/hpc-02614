from experiment import RunConfig, experiment

config = RunConfig(gauss_seidel=True, n=100, iter_max=10_000, tolerance=0)

experiment([config], "data.csv")
