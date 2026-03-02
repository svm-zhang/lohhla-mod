<h1>
    lohhlamod
</h1>


## Introduction

`lohhlamod` is the original [LOHHLA](https://doi.org/10.1016/j.cell.2017.10.001) algorithm re-engineered in modern style that

- offers additional features/metrics for better interpretation
- offers runtime speedup
  - separates HLA realignment from LOH detection
  - maximizes vectorized operations whenever possible
  - uses data.table for efficient data processing
- provides better output layout and makes intermediate result per allele in RDS format accessible
- provides proper packaging for ease use in multi-user HPC environment
- removes hardcoded path presets for better code maintenance

## Installation

It is recommended to use docker for building and running `lohhlamod`.

```
docker compose build --pull
```

## Run LOH analysis

Please refer to the [documentation](https://svm-zhang.github.io/lohhla-mod/) for details.

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

