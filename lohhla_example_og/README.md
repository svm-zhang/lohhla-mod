To run `lohhlamod` on the original example data:

```
hlalohReforged --tbam example_tumor.hla.realn.ready.bam \
    --nbam example_normal.hla.realn.ready.bam \
    --hlaref example.patient.hlaFasta.fa \
    --tstates tstates.tsv \
    --outdir $PWD/example \
    --subject example
```

Note that the input alignment files used here are not the original BAM files provided by the `lohhla` package. I use the following steps:

1. HLA typing on the `example_BS_GL_sorted.bam` using `polysolvermod`
2. Realign the `example_tumor_sorted.bam` on the typing results using `polysolvermod`

The resulting `*.ready.bam` files for both normal and tumor samples are used as input to `lohhlamod` command shown here.

The `example.loh.res.tsv` gives you the LOH result with this original example data.