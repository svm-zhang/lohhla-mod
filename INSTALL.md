`lohhla-mod` uses [mamba](https://github.com/mamba-org/mamba) and [boa](https://github.com/mamba-org/boa) as package manager and builder.

## Linux (linux-64)

To install `lohhla-mod` on a LINUX platform, simply run the following:

```
cd "$lohhlamod_repo"
boa build .     # this builds polysolvermod as a local tarball
mamba create -n lohhla
mamba install -n lohhla --use-local lohhlamod
```

Upon successful installation, you should be prompted with help message like below:

```
usage: lohhlamod
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
  --threads INT    Specify the number of threads (16)
```

## MacOS (osx-arm64, bash)

Installation on OSX-64 requires a few additional steps prior to running `boa build`:

- You need to switch to `bash` shell
- Set `subdirs` to `osx-64` in your `.condarc` or `.mambarc` config file

## Manual

It is also possible to install `lohhla-mod` without any package manager and builder. The `recipe.yaml` file defines run-time dependencies. You can install them on your local machine/environment as you see fit. Once everything installed, running the code below should get you going:

```
bin_dir="${lohhlamod_repo}/bin"
mkdir "$bin_dir"
cp "${lohhlamod_repo}/R/lohhlamod.R" "${bin_dir}/lohhlamod"
cp "${lohhlamod_repo}/R/lohhlaplot.R" "${bin_dir}/lohhlaplot"
export R_LIBS="$lohhlamod_repo:$R_LIBS"
export PATH="$bin_dir:$PATH"
```

If you want to run `lohhla-mod` everytime you open a session, simply add the last two lines to your `~/.bash_profile` file.
