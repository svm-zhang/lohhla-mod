
We use [mamba](https://github.com/mamba-org/mamba) as the package manager, and [boa](https://github.com/mamba-org/boa) to build it.

To install mamba, please follow the [instruction](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-mac-os--linux) from [Miniforge3](https://github.com/conda-forge/miniforge):

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh
```

Start a new shell session after the installation finishes. `mamba` binary now should be accessible and you should be able to run:

```
mamba -h
```

Before we move onto the next step, it is recommended to set up the `CONDA_BLD_PATH` environment variable, and it is preferrable to set it in your `.bash_profile` like below. This gives `boa` a default location to output packages (see `--output-folder` option in `boa build` command). 

```
echo "export CONDA_BLD_PATH=${MINIFORGE3_PREFIX}/conda-bld" >> ~/.bash_profile
```

Now simply install `boa`:

```
mamba install boa -c conda-forge
```

if you do not have `git` installed at this moment,

```
mamba install git
```

## Linux (linux-64)

On Linux machine, first clone the hlalohReforged repo:

```
git clone https://github.com/svm-zhang/hlalohReforged.git
```

Then, move to the local directory where you just clone the package and run `boa` to build the package:

```
cd hlalohReforged

boa build . 
```

`boa`, by default, should output the package tarball under `${MINIFORGE3_PREFIX}/conda-bld/linux-64`. `${MINIFORGE3_PREFIX}` is where Miniforge3 is installed, and `linux-64` is the target-platform. `linux-64` should be the most common architecture used in the bioinfo field. __Other linux architecture platform is not supported__.

Next, build a conda environment and install the `hlalohreforged` package you just build using mamba:

```
mamba create -n hlalohreforged      # create an environment

mamba install -n hlalohreforged --use-local hlalohreforged      # hlalohreforged is a local package

```

Last, activate the environment with the `hlalohReforged` installed and start running:

```
mamba activate hlalohreforged

hlalohReforged -h
```

Upon successful installation, you should be prompted with help message like below:

```
usage: /home/simo/opt/miniforge3/envs/lohreforged/bin/hlalohReforged
       [-h] --subject STR --tbam FILE --nbam FILE --hlaref FILE
       [--tstates FILE] --outdir DIR [--min_cov INT] [--min_necnt INT]
       [--threads INT]

options:
  -h, --help       show this help message and exit
  --subject STR    Specify the subject ID
  --tbam FILE      Specify the tumor bam file
  --nbam FILE      Specify the normal bam file
  --hlaref FILE    Specify HLA reference sequence
  --tstates FILE   Specify file includeing tumor purity and ploidy
  --outdir DIR     Specify the output directory
  --min_cov INT    Specify the minimum coverage at mismatch sites (30)
  --min_necnt INT  Specify the minimum number of diff events allowed for reads
                   mapping to HLA alleles (1)
```

## MacOS (osx-arm64, bash)

The general processing of building `hlaohreforged` on MacOS is similar as above. The only difference is that we need to tell `boa` to use `osx-64` as the target platform, rather than `osx-arm64`. This is because `bioconda` does not support `osx-arm64` yet. First, run the following command to make `osx-64` subdirectory.

```
conda config --env --set subdir osx-64
```

Then, we specify target platform in `boa build` command.

```
boa build --target-platform osx-64 . 

```

All other steps follow the same procedure as decribed in the Linux section above.

## Manual

It is possible to install `hlalohReforged` without any package manager and builder. The `recipe.yaml` file defines run-time dependencies. You can install them on your local machine/environment as you see fit. Once everything installed, you need to set one more environment variable to run `hlalohReforge`:

```
export R_LIBS=${Path to the parent directory of hlalohReforged repo}:$R_LIBS
```

You can also set this in your `~/.bash_profile` file. 
