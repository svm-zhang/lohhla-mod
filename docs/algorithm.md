# Algorithm & interpretation

## Key Enhancements over the Original LOHHLA

While the original `LOHHLA` provides a robust foundation for LOH detection,
`lohhlamod` introduces several architectural and statistical refinements
designed to reduce false positives and improve robustness in real-world
clinical datasets.

### Multi-metric evaluation

The original algorithm primarily relies on a copy number (CN) threshold (<0.5)
and a unique P-value (`PVal_unique` < 0.01). To address observed overcalling—particularly
in samples with noisy coverage—`lohhlamod` introduces additional bin-level statistics:

- `Pct_CN_Diff_Supporting_Bins`: Measures the proportion of 150bp bins that
support a significant CN difference between the two alleles.
A true LOH event typically spans the entire length of the allele. A suggested
starting threshold for high-confidence calls is 75%.

    ???+ Note "Threshold value"
        In real-world datasets, it is not always possible to observe a 100% bins
        supporting a copy number difference due to the technical difficulties
        imposed by the highly-polymorphic nature of MHC regions. Utilizing
        this value can help filter out suspicious false positive LOH calls.

- **Allele-Specific Loss Metrics** (`Pct_A1/A2_Loss_Supporting_Bins`): These metrics
use a one-sample t-test of allelic logR (H0: μ=−1) to determine if an allele shows a
statistically significant loss. This is particularly useful for
cases where estimated CN for both alleles drop below 0.5 or even become negative.

    ???+ Note "Why use direct test"
        Relying solely on copy number estimates can be risky with noisy coverage
        data. Having direct statistical tests for allelic loss provides a more
        grounded way to identify true deletions.

    
### BAF corrected for allelic capture bias

In many assays, one HLA allele may be captured more efficiently than its
counterpart due to probe affinity or sequence context. `lohhlamod` estimates
this "capture bias" from the **normal sample** and applies a correction factor
to the observed BAF in the tumor.

???+ Note "Plausibility check"
    If the correction results in a BAF >1.0 (which can happen if the capture
    bias and tumor observation directions differ), lohhlamod forces the BAF
    to a value of 1.0 to maintain biological plausibility. When this happens,
    it is recommended to:

    - Verify the HLA typing result.
    - Check for potential sample swapping.

### Global depth corrector

Log-ratio (LogR) calculation depends on normalizing depth differences between
tumor and normal samples. In an ideal world, a neutral diploid state
should results in a LogR value of 0.

In reality, however, the diploid baseline often shifts away from zero due to
technical challenges. Furthermore, a global, library-wise depth corrector
may not reflect the specific normalization needed locally across the HLA genes.

Unlike the original tool which calculates a global corrector across the entire
BAM file, `lohhlamod` estimates the corrector locally from the HLA genes. By
applying the high-quality alignment restrictions (via `--min-necnt`) to the
depth corrector as it does to the coverage calculation, the tool ensures the
normalization baseline is more representative of the local observation.

## Suggested Interpretation Workflow

When reviewing results, I recommend following this hierarchy:

1. Check Bin Consistency
    Look at `Pct_CN_Diff_Supporting_Bins`. If the
    proportion is low, the imbalance is likely localized noise rather than
    a true LOH event.

    ???+ Note "Be careful"
        This metric identifies the percentage of bins supporting a copy number difference.
        However, an amplification event can also produce a high value. Therefore, it is
        recommended not to use this metric in isolation.

2. Evaluate Specific Loss
    Check `Pct_A1_Loss_Supporting_Bins` and `Pct_A2_Loss_Supporting_Bins`.

    - **In the case of LOH**: Either metric is expected to be high (see the [HLA-A result](https://github.com/svm-zhang/lohhla-mod/blob/main/simulation/s6/output/s6.loh.res.tsv)
        in simulated `s6` case).

    - **In the case of Homozygous Loss**: Both metrics are expected to be high (see the [HLA-C result](https://github.com/svm-zhang/lohhla-mod/blob/main/simulation/s7/output/s7.loh.res.tsv)
        in simulated `s7` case).

3. Review LogR and BAF
    Ensure the median LogR estimates and BAF shift align
    with the bin-level statistics. A true loss should show a consistent shift
    away from **subject-specific neutral states** across both metrics.

4. Validate CN Estimates
    - If the estimated copy number for both alleles are negative or near zero,
        review your coverage data via inspecting tables saved in the `.rds` files
        and/or checking the diagnostic plots.
    - If the estimated copy number falls outside the calculated lower and upper
        bounds, it suggests the data may be skewed and the statistical assumptions
        may be violated.

    ???+ Note "Can you trust copy number estimates"
        Copy number estimates are sensitive to noise in the coverage data.
        While the estimation is reliable in well-profiled assays (as
        demonstrated in the "best-case" simulation scenarios),
        `lohhlamod` is currently best used as a qualitative measure to determine
        if a LOH event has occurred.

5. Visual Confirmation: Always inspect the diagnostic plots generated by
    `lohhlaplot` command, especially the [paired allelic coverage distribution plot](./visual/#paired-allelic-coverage-distribution).
    Visual evidence remains a reliable validator.

## Hidden cutoffs

The following parameters are currently hardcoded but may be exposed to the CLI in future releases:


| Parameter | Value | Description |
| -------- | ----- | ----- |
|  Bin size |  150bp  | The sliding window size for coverage calculation  |
|  Min. mismatches  |  5  | Minimum number of mismatches required between alleles to attempt LOH detection. |
|  Alpha  |  0.01  | Significance threshold for statistical tests.  |
|  Gamma  |  1  | Smoothing parameter for copy number estimation.  |
|  Base quality  |  20  | Minimum base quality required for `samtools mpileup`.  |
|  SAM flag  |  3584 / 2  | Excludes **UNMAP**, **QCFAIL**, **DUP** (3584); includes **PROPER_PAIR** (2).  |



