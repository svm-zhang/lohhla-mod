# lohhlamod


## Introduction

`lohhlamod` is the original
[LOHHLA](https://doi.org/10.1016/j.cell.2017.10.001) algorithm re-engineered in
modern style that

- offers additional features/metrics for better interpretation
- offers runtime speedup
  - separates HLA realignment from LOH detection
  - maximizes vectorized operations whenever possible
- provides better output layout and makes intermediate result per allele in RDS
  format accessible
- provides proper packaging for ease use in multi-user HPC environment
- removes hardcoded path presets for better code maintenance

## Installation

It is recommended to use docker for building and running `lohhlamod`.

=== "docker"
```
docker compose build --pull
```

Add `--no-cache` option for a refresh build.

???+ "Verify installation"
    After installation, I recommend first running `lohhlamod` using the sample
    simulation data provided in the repository, and verify
    that your output matches the expected results. The expected results can be
    found under the `output/` subdirectory within each simulated case. Refer to
    the [Dataset](./dataset.md) for details about simulation.


## Citation

Please cite the original [HLALOH](https://doi.org/10.1016/j.cell.2017.10.001)
paper, its
[Bitbucket](https://bitbucket.org/mcgranahanlab/lohhla/src/master/) repository, and
the [lohhlamod](https://github.com/svm-zhang/lohhla-mod) repository.


## License

- `lohhlamod` fully respects all [LICENSE
  requirments](https://bitbucket.org/mcgranahanlab/lohhla/src/master/) imposed
  by the original `LOHHLA` tool.
- `lohhlamod` is free to use for all non-commercial parties.

## Disclaimer

- I, by no means, intent to overtake the origianl idea, implementation, and
  copyright of the original `LOHHLA` algorithm.
- This repo does not distribute `LOHHLA` package, as well as all its dependent
  third-party parties that are under commercial licenses.
- `lohhlamod` does not necessarily produce identical result as `LOHHLA`,
  including but not limited to estimates of copy numbers, logR, BAF, plots,
  etc.
- Please interpret result at your own discretion when using `lohhlamod`.
  I simulated a few datasets to demonstrate the performance of `lohhlamod` in
  the most ideal settings.

