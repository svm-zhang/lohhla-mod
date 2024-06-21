## Case: Subject s6

| HLA Gene | A1_CN | A2_CN |
| -------- | ----- | ----- |
|  HLA-A   |   1   |   0   |
|  HLA-B   |   1   |   1   |
|  HLA-C   |   3   |   1   |

LOH event occurrs at HLA-A gene. In addition, subject `s6` also has an amplification of one of the alleles at HLA-C gene.

Please use this [file](./s6/s6_polysolvermod/s6_n.hla.fasta.fai) for details about HLA genotypes of subject `s6`.

```
hlalohReforged --tbam s6_t.hla.realn.ready.bam \
    --nbam s6_n.hla.realn.ready.bam \
    --hlaref s6_n.hla.fasta \
    --tstates tstates.tsv \
    --outdir $PWD/s6/s6_lohhlamod \
    --subject s6
```

The command above generates `s6.loh.res.tsv` LOH result within the output directory you specified.

You can apply the same command to all other cases below, with different BAM and HLA reference files.

## Case: Subject s7

| HLA Gene | A1_CN | A2_CN |
| -------- | ----- | ----- |
|  HLA-A   |   1   |   1   |
|  HLA-B   |   2   |   1   |
|  HLA-C   |   0   |   0   |

Subject `s7` loses both alleles of HLA-C gene.

## Case: Subject s1

| HLA Gene | A1_CN | A2_CN |
| -------- | ----- | ----- |
|  HLA-A   |   2   |   0   |

This is a scenario of a copy neutral LOH event at HLA-A gene. Note that this and the next are the most simplified cases I first simulated for developing `lohhlamod`.

## Case: Subject s2, a non-LOH case

| HLA Gene | A1_CN | A2_CN |
| -------- | ----- | ----- |
|  HLA-B   |   1   |   1   |

## How to simulate

### Simulate normal sample
1. Create a _in silico_ normal individual with HLA alleles at HLA gene loci
2. Make a HLA reference file with the selected alleles in Fasta format
3. Simulate sequencing reads from the HLA reference using [dwgsim](https://github.com/nh13/DWGSIM/tree/main). The example command below generates 10k paired-end reads of 150bp from the HLA reference. Please refer to `dwgsim` documentation for details about all the command line options.
```
dwgsim -d 300 -s 20 -N 10000 -1 150 -2 150 -R 0.02 -S 2 -H -o 1 -e 0.00109 -E 0.00109 "$HLA_REF" "$OUT_PREFIX"
```
4. Align simulated reads against human genome, sort and index the resulting BAM file
5. Run `polysolvermod` to genotype HLA alleles using the BAM

### Simulate tumor sample
Now switch gear to simulate tumor data from the normal data.
1. Decide copy number of each allele at each HLA gene locus
2. Make a modified HLA reference from the normal one by honoring the copy numbers. Note that if you make a LOH event, the deleted alleles should not in the reference, instead you need to create a separate Fasta file for the delete alleles
3. Simulate sequencing reads from the modified tumor HLA reference file using `dwgsim` command similar as above. The example command below generates 20k paired-end reads of 150bp with higher error rates and shorter fragment templates (trying to mimic tumor sample).
```
dwgsim -d 280 -s 20 -N 20000 -1 150 -2 150 -R 0.02 -S 2 -H -o 1 -e 0.002 -E 0.002 "$Tumor_HLA_REF" "$OUT_PREFIX"
```
4. Simulate sequencing reads for the delete alleles separately at lower coverage to mimic the loss event. Adjust the `-N` option in the `dwgsim` command can achieve this. 
5. Combine all simulated sequencing reads and do step 4 as in the normal case
6. Run `polysolvermod` with the `--realn_only` option to get the BAM file for `lohhlamod`. The HLA reference is the output reference from step 5 in the normal case

