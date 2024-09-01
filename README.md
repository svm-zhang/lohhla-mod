<h1>
    lohhla-mod
</h1>

<!-- TOC start (generated with https://github.com/derlin/bitdowntoc) -->

- [Introduction](#introduction)
- [Installation](#installation)
- [Command line](#command-line)
- [Prepare Input](#prepare-input)
  - [BAM](#bam)
  - [HLA reference](#hla-reference)
  - [Estimated tumor ploidy and purity (--tstates)](#estimated-tumor-ploidy-and-purity---tstates)
- [Explain Output](#explain-output)
- [Visualize Coverage, LogR, and BAF](#visualize-coverage-logr-and-baf)
- [Simulation data](#simulation-data)
- [Original lohhla test data](#original-lohhla-test-data)
- [Key differences from the OG lohhla algorithm](#key-differences-from-the-og-lohhla-algorithm)
  - [Additional metrics for better interpretation](#additional-metrics-for-better-interpretation)
  - [BAF corrected for allelic capture bias](#baf-corrected-for-allelic-capture-bias)
  - [Global depth corrector](#global-depth-corrector)
- [Suggested Interpretation using lohhla-mod](#suggested-interpretation-using-lohhla-mod)
- [Hidden cutoffs](#hidden-cutoffs)
- [License](#license)
- [Disclaimer](#disclaimer)
- [Citation](#citation)

<!-- TOC end -->

## Introduction

`lohhla-mod` is the original [LOHHLA](https://doi.org/10.1016/j.cell.2017.10.001) algorithm re-engineered in modern style that

- offers additional features/metrics for better interpretation
- offers runtime speedup
  - separates HLA realignment from LOH detection
  - maximizes vectorized operations whenever possible
  - uses data.table for efficient data processing
- provides better output layout and makes intermediate result per allele in RDS format accessible
- provides proper packaging for ease use in multi-user HPC environment
- removes hardcoded path presets for better code maintenance

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

`lohhla-mod` does not do realignment like the original `LOHHLA` program. To get the required BAMs, you can use [polysolverMod](https://github.com/svm-zhang/polysolverMod), a re-engineered HLA typing tool based on `Polysolver`.

To get BAM for normal sample, you can simply follow the [example](https://github.com/svm-zhang/polysolverMod?tab=readme-ov-file#quick-start) with your data. `polysolverMod` generates the realigned BAM with suffix `ready.bam`, that is ready for detecting LOH.

Please follow this [guide](https://github.com/svm-zhang/polysolverMod?tab=readme-ov-file#scenario-detecting-loh-from-paired-tumor-and-normal-samples) specifically for getting realigned BAM for tumor sample.

### HLA reference

This refers to the sample-level HLA reference with specific typed alleles for the sample. The reference file is available after you successfully run `polysolverMod` on the normal sample.

### Estimated tumor ploidy and purity (--tstates)

`lohhla-mod` uses estimated ploidy and purity for inferring allelic copy number. Ploidy and purity estimates can be obtained from many CNV algorithms. In cases where there is no paired normal sample or reference panel available such that you cannot get the estimates, `lohhla-mod` allows this option to be optional by using a default value of `ploidy = 2` and `purity = 0.5`. The default values have not made customizable from command line at the moment.

An example of the file provided to `--tstates` looks like below:

| SampleID | TumorPloidy | TumorPurityNGS |
| -------- | ----------- | -------------- |
| s1_t     | 2.33        | 1              |

## Explain Output

`lohhla-mod` dumps all results under output specified by `--outdir`.

- `*.filt.bam`: filtered alignment result by minimum allowed mismatch events specified by `--min_ecnt` option. The filtered BAM files are used for the final LOH detection
- `$subject.loh.res.tsv`: LOH main result, each row per allele (see below for column schema)
- `*.rds`: serialized file with intermediate data tables, one per allele. These `rds` files are used for getting plots by running `lohhlaplot` command.

The columns in the LOH result are defined as follows:

- `HLA_A1_CN`: estimated copy number (CN) for allele 1, the upper and lower estimates are named with `Upper` and `Lower` suffix
- `HLA_A2_CN`: estimated CN for allele 2
- `Pct_CN_Diff_Supporting_Bins`: percentage of bins supporting a significant CN difference b/w allele 1 and 2
- `HLA_A1_Median_LogR`: median estimates of log-ratio of tumor versus normal for allele 1
- `HLA_A2_Median_LogR`: same as above but for allele 2
- `HLA_A1_MM_Median_LogR`: median estiamtes of log-ratio of tumor versus normal at mismatch sites for allele 1
- `HLA_A2_MM_Median_LogR`: same as above but for allele 2
- `MM_LogR_Paired_Pvalue`: p-value from paired test for log-ratio difference b/w 2 alleles at mismatch sites
- `Median_BAF`: median estimate of b-allele frequency across mismatch sites
- `Num_MM`: number of mismatches b/w 2 alleles
- `Num_Bins`: total number of bins of size 150bp across the pairwise alignment b/w 2 alleles
- `Num_MM_Bins`: number of bins with mismatch sites
- `Pct_A1_Loss_Supporting_Bins`: percentage of bins supporting a significant CN loss for allele 1
- `Pct_A2_Loss_Supporting_Bins`: same as above but for allele 2
- `HLAGene`: HLA gene locus
- `HLA_A1`: allele 1
- `HLA_A2`: allele 2

## Visualize Coverage, LogR, and BAF

`lohhla-mod` provides a separate command `lohhlaplot` to generate a set of plots per HLA gene. `lohhlaplot` uses the `rds` files as backend data for visualization. To get the plots provided by this package, simply run:

```
lohhlaplot --sample "$subject" \
    --loh_res "$loh_res_file" \
    --loh_dir "$loh_outdir"
```

The command above creates a sub-folder within the `--loh_dir` folder. Each HLA gene gets a set of plots:

1. coverage distribution across the entire length of allele 1 and 2 in the normal sample
   ![s6.hla_a.n_dp.png](./simulation/s6/s6_lohhlamod/s6_plots/hla_a.n_dp.png)
2. coverage distribution across the entire length of allele 1 and 2 in the tumor sample
   ![s6.hla_a.t_dp.png](./simulation/s6/s6_lohhlamod/s6_plots/hla_a.t_dp.png)
3. paired tumor and normal coverage distribution in log10 scale across the entire length of allele 1 and 2. Tumor coverage in this plot is corrected for tumor and normal depth difference. This plot should inform you whether or not there is a LOH event occurring in one of the alleles.
   ![s6.hla_a.tn_dp.png](./simulation/s6/s6_lohhlamod/s6_plots/hla_a.tn_dp.png)
4. distribution of log-ratio of tumor versus normal across entire length of allele 1 and 2. The dashed line represents median estimate of logR for each allele
   ![s6.hla_a.logR.png](./simulation/s6/s6_lohhlamod/s6_plots/hla_a.logR.png)
5. distribution of BAF
   ![s6.hla_a.baf.png](./simulation/s6/s6_lohhlamod/s6_plots/hla_a.baf.png)

## Simulation data

I provided a few _in silico_ datasets to mimic certain scenarios. But it is important to note that these are idealized datasets that provide guidance only. Please refer to the [README](./simulation/README.md) for more details.

## Original lohhla test data

The original example data from the `LOHHLA` algorithm is also available in this repo. I used `polysolverMod` to get the BAM files and HLA reference. Then `lohhla-mod` was used to generate the LOH result. You can use this data to test `lohhla-mod` prior to running on your own data. Note that there is LOH result for HLA-A gene as it is the only HLA gene provided in the example.

## Key differences from the OG lohhla algorithm

### Additional metrics for better interpretation

The original `lohhla` calls a LOH event when following conditions are met:

- A copy number < 0.5 (Step 5 section under Method in the paper)
- `PVal_unique < 0.01` (Step 5 section under Method in the paper)
- The lost allele is the one with the lower median logR at mismatch sites (line 1404-1405 in `LOHHLAscript.R`)

These are all valid and good choices to have good LOH detection accuracy. However, I found these parameters can still lead to overcalling of LOH events, at least on the datasets I was working on in the past.

`lohhla-mod` offers a few more metrics to help (hopefully):

- `Pct_CN_Diff_Supporting_Bins`: this metric tells you how many bins support a significant CN difference b/w the 2 alleles. If there is a LOH event (e.g. 2/0 or 1/0), it is expected that the event spans across a large portion of the allele length. Personally, I find `75%` is a starting point to tune this value
- `Pct_A1_Loss_Supporting_Bins` and `Pct_A1_Loss_Supporting_Bins`: these 2 metrics tell you how many bins support a CN loss for A1 and A2 alleles, respectively. By "CN loss", it is coded as a one-sample t-test of allelic logR against `mu=-1`. The motivation behind having these 2 metrics is to handle cases (I encountered a lot) where estimated CNs for both alleles are less than `0.5`, or are even negative. Personally, I am not knowledgeable enough to have an approximate number on the proprotion of individuals in a clinical study who lose both copies of a HLA gene. But I was handed LOH results from running `lohhla` that had both CN estimates less than 0.5 or 0
- `MM_LogR_Paired_Pvalue`: this metric is the same as the `PVal` (legacy) metric, rather than `PVal_unique`. I do not quite get the gist of the `unique` concept underlying `PVal_unique`. And quite honest, I personally do not think it will greatly help avoid to overcall LOH. That is why I choose to use a simple solution

### BAF corrected for allelic capture bias

When one allele has a higher frequency than its counterpart at a site, it can also be due to the fact that the assay (weblab) captures one allele better than the other. `lohhla-mod` estimates capture bias from the normal sample, and corrects for the observed BAF in the tumor.

Note that the correction can lead to a BAF larger than 1.0. This happens when at some sites, the capture bias corrector does not have the same direction as the observation in tumor. `lohhla-mod` forces BAF to be a value of 1 in such cases.

### Global depth corrector

logR is calculated by correcting for the depth difference b/w tumor and normal libraries. In an idealized world, this correction should always make logR a value of 0 in a diploid state. In reality, however, I have seen considerable number of cases where median logR are negative for both alleles, which in turn leads to underestimation of copy number. In such cases, negative logR are not necessarily indicative of copy number loss. Rather, it means i) the diploid state of logR is shift from zero; ii) the global depth corrector might not reflect the observation across HLA genes.

`lohhla-mod` tries to make the global estimator reflecting what happens locally. And because the local coverage calculation uses the `--min_ecnt` option, `lohhla-mod` simply adds the restriction when calculating the global corrector. The solution helps in certain cases, but does not work for everything as much as I would like.

## Suggested Interpretation using lohhla-mod

1. First look at `Pct_CN_Diff_Supporting_Bins`, if you do not see a high proportion of bins supporting a CN difference, there is likely not a LOH event. Note that an amplification event can also have a high number for this metric
2. Next look at `Pct_A1_Loss_Supporting_Bins` and `Pct_A1_Loss_Supporting_Bins` metrics. When LOH happens, either A1 or A2 should have a high number (check out the HLA-A case in simulated `s6` case). When both alleles are lost, both metrics should be high (check out the HLA-C case in simulated `s7` case)
3. Then look at median logR estiamtes and median BAF metrics
4. Last look at copy number estimates for both alleles. Note that if you see the estimated CN is outside of lower and upper bounds, it means you have a skewed data and the assumption of t test used for the estimate is probably being violated. This in turn also means that you need to take a closer look at your data
5. I would also suggest check out the plot, especially the `tn_dp` one.

## Hidden cutoffs

There are a few pre-defined and non-customizable cutoffs used in `lohhla-mod`. These cutoffs can be exposed to command line if needed in the future:

- bin size: 150bp
- minimum number of mismatches b/w 2 alleles required: 5 (below which no LOH detection will be attempted)
- p value: 0.01
- gamma: 1
- parameters used for running samtools mpileup
  - minimum base quality: 20
  - exclude flag: 3584 [UNMAP, QCFAIL, DUP]
  - include flag: 2

## License

- `lohhla-mod` fully respects all [LICENSE requirments](https://bitbucket.org/mcgranahanlab/lohhla/src/master/) imposed by the original `LOHHLA` tool.
- I am currently checking with the original authors of `LOHHLA` package, and will update license file once I get confirmation. For now, `lohhla-mod` is free to use for all non-commercial parties.

## Disclaimer

- I, by no means, intent to overtake the origianl idea, implementation, and copyright of the original `LOHHLA` algorithm.
- This repo does not distribute `LOHHLA` package, as well as all its dependent third-party parties that are under commercial licenses.
- `lohhla-mod` does not necessarily produce identical result as `LOHHLA`, including but not limited to estimates of copy numbers, logR, BAF, plots, etc.
- Please interpret result at your own discretion when using `lohhla-mod`. I simulated a few datasets to demonstrate the performance of `lohhla-mod` in the most ideal settings.

## Citation

Please cite the original [HLALOH](https://doi.org/10.1016/j.cell.2017.10.001) paper and its [Bitbucket](https://bitbucket.org/mcgranahanlab/lohhla/src/master/) repository.

If you use `lohhla-mod`, please kindly cite this github repo as well.
