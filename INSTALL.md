
We use [mamba](https://github.com/mamba-org/mamba) as the package manager, and [boa](https://github.com/mamba-org/boa) to build it.

To install mamba, please follow the [instruction](https://github.com/conda-forge/miniforge?tab=readme-ov-file#unix-like-platforms-mac-os--linux) from [Miniforge3](https://github.com/conda-forge/miniforge):

```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh)
```

Once the installation finishes, start a new shell session and `mamba` binary should be accessible. Next, install `boa` by running:

```
mamba install boa -c conda-forge
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

`boa`, by default, should output the package tarball under `${MINIFORGE3_PREFIX}/conda-bld/linux-64`. `${MINIFORGE3_PREFIX}` is where Miniforge3 is installed, and `linux-64` is the target-platform.

If it is not there, please check the environment variable `$CONDA_BLD_PATH`. To set it to where the installation location of Miniforge3, simply do:

```
echo "export CONDA_BLD_PATH=${MINIFORGE3_PREFIX}/conda-bld" >> ~/.bash_profile
```

Next, build a conda environment to run `hlalohReforged` using mamba:

```
mamba create -n hlalohreforged      # create an environment

mamba install -n hlalohreforged --use-local hlalohreforged      # hlalohreforged is a local package

```

Last, activate the environment with the `hlalohReforged` installed and start running:

```
mamba activate hlalohreforged

hlalohReforged -h
```

The last command should print out the help message when everything works as expected.

## MacOS (osx-arm64)

The general processing of building `hlaohreforged` on MacOS is similar as above. The only difference is that we need to tell `boa` to use `osx-64` as the target platform, rather than `osx-arm64`. This is because `bioconda` does not support `osx-arm64` yet. First, run the following command to make `osx-64` subdirectory.

```
conda config --env --set subdir osx-64
```

Then, we specify target platform in `boa build` command.

```
boa build --target-platform osx-64 . 

```

All other steps follow the same procedure as decribed in the Linux section above.


