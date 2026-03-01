# Getting started

## Quick example

Using the simulation data, [s6](https://github.com/svm-zhang/lohhla-mod/tree/main/simulation), provided in the repository as an example:

``` bash
docker compose run --rm lohhlamod \
    lohhlamod --tbam ./simulation/s6/input/s6_t.hla.realn.ready.bam \
    --nbam ./simulation/s6/input/s6_n.hla.realn.ready.bam \
    --subject s6 \
    --hlaref ./simulation/s6/input/s6_n.hla.fasta \
    --tstates ./simulation/s6/input/tstates.tsv \
    --outdir ./s6_test
```

The output folder `s6_test` contains the following results. Details descriptions
of each file can be found on the [Explain Output](./output.md) page.

``` { .sh }
s6_test/
├─ hla_a.rds
├─ hla_b.rds
├─ hla_c.rds
├─ s6.loh.res.tsv
├─ s6_n.filt.bam
├─ s6_n.filt.bam.bai
├─ s6_t.filt.bam
├─ s6_t.filt.bam.bai
```


## Input preparation

`lohhlamod` is designed as a specialized, high-performance module focused
exclusively on the detection of HLA Loss of Heterozygosity (LOH).

While the original LOHHLA framework included an internal HLA typing and
realignment routine—which users could optionally skip—lohhlamod explicitly
removes these components to provide a leaner, more modular footprint. This
architectural choice treats HLA typing as an upstream prerequisite rather
than an internal step, reflecting the reality that most modern bioinformatics
pipelines already have a preferred HLA typing method in place.

### Required inputs

To run `lohhlamod`, you must provide:

- HLA-aligned BAMs: Alignments against subject's specific HLA reference
    for both normal and tumor samples.
- HLA reference: The subject-specific HLA alleles in FASTA format.
- Tumor states: A TSV file containing tumor ploidy and purity estimates.

### HLA reference and alignment

[mhcflow](https://github.com/svm-zhang/mhcflow), is a re-engineered HLA typing
tool based on `Polysolver`, which is optimized to produce the specific
alignments required for downstream LOH analysis. You can find detailed preparation
instructions [here](https://svm-zhang.github.io/mhcflow/loh/).

???+ Tip "Input compatibility"
    `lohhlamod` expects BAM files where reads have been re-aligned to the
    patient's inferred HLA alleles (e.g., HLA-A, B, and C). If you are using a
    custom pipeline instead of `mhcflow`, ensure your HLA
    reference contains only the subject-specific alleles and that you align
    reads against this specific reference for both tumor and normal samples.


### Estimated tumor ploidy and purity (--tstates)

`lohhlamod` uses estimated ploidy and purity to infer allelic copy
number. These estimates can be obtained from most standard CNV algorithms.

???+ Note "Default ploidy and purity values"
    In cases where estimates are unavailable, lohhlamod will fall back to default
    values of `ploidy = 2` and `purity` = 0.5. Currently, these defaults
    are hardcoded and cannot be modified via the command line.

???+ Tip "Format"
    The file passed to --tstates must be a tab-delimited file. You can refer
    to the `tstates.tsv` files provided under `simulation/`
    directory in this repository for exact formatting details.

Example:


| SampleID | TumorPloidy | TumorPurityNGS |
| -------- | ----------- | -------------- |
| s1_t     | 2.33        | 1              |


## Command line interface

```text
usage: lohhlamod
       [-h] --subject STR --tbam FILE --nbam FILE --hlaref FILE
       [--tstates FILE] --outdir DIR [--min-cov INT] [--min-necnt INT]
       [--threads INT]

options:
  -h, --help       show this help message and exit
  --subject STR    Specify the subject ID
  --tbam FILE      Specify the tumor bam file
  --nbam FILE      Specify the normal bam file
  --hlaref FILE    Specify HLA reference sequence
  --tstates FILE   Specify file includeing tumor purity and ploidy
  --outdir DIR     Specify the output directory
  --min-cov INT    Specify the minimum coverage at mismatch sites (30)
  --min-necnt INT  Specify the minimum number of diff events allowed for reads
                   mapping to HLA alleles (1)
  --threads INT    Specify the number of threads (16)
```

