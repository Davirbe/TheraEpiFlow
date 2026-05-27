"""TheraEpiFlow validation suite.

Benchmark + reproducibility experiments for the paper. Runs the production
pipeline against disposable projects, captures per-step timing and per-stage
epitope counts, and renders publication-ready figures.

Strictly serial — never introduces additional parallelism so wall-clock
measurements stay valid. The pipeline's internal NetMHCpan + MHCFlurry
ThreadPoolExecutor is preserved as production-equivalent.
"""
