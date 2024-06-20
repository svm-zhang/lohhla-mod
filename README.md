# lohhla Modern

`lohhlamod` is the original [lohhla](https://doi.org/10.1016/j.cell.2017.10.001) HLA loss of heterozygosity detection algorithm re-engineered in modern style.  

## Installation

Please refer to [INSATLL](INSTALL.md) for details.

## Command line

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

## Prepare Input

### BAM

`lohhlamod` does not do realignment like the original `lohhla` program. To get the required BAMs, you can use [polysolverMod](https://github.com/svm-zhang/polysolverMod), another re-engineered HLA typing tool based on `polysolver`.

To get BAM for normal sample, you can simply follow the [example](https://github.com/svm-zhang/polysolverMod?tab=readme-ov-file#quick-start) and swap with your data. `polysolverMod` generates the realigned BAM with suffix `ready.bam`, that is ready for detecting LOH.

Please follow this [guide](https://github.com/svm-zhang/polysolverMod?tab=readme-ov-file#scenario-detecting-loh-from-paired-tumor-and-normal-samples) specifically for getting realigned-BAM for tumor sample. 

## Explain Output


## Key differences from the OG lohhla algorithm


## Disclaimer

I, by no means, try to overtake the origianl idea and implementation of `polysolver` algorithm. This repo opens to all non-commercial researchers and projects. My mere purpose is to make `polysolver` better, faster, more versatile, if not more accurate.

## Citation

Please cite the original [polysolver](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4747795/) paper.

If you use `polysolvermod`, please cite this github repo as well.