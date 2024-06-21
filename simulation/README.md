## Case: Subject s6

| HLA Gene | A1_CN | A2_CN |
| -------- | ----- | ----- |
|  HLA-A   |   1   |   0   |
|  HLA-A   |   1   |   1   |
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

## Case: copy-neutral LOH

## Case: Non-LOH


## Case: 1/1;2/1;0/0

Simulate a case where HLA-C gene locus is entirely deleted
```
hlalohReforged --tbam s7_t.hla.realn.ready.bam --nbam s7_n.hla.realn.ready.bam s7_n.hla.fasta --tstates tstates.tsv --outdir $PWD --subject s7 --min_necnt 2 --corrector global
```


## How to simulate

I use [dwgsim](https://github.com/nh13/DWGSIM/tree/main) to simulate a few datasets as guidance set up expectation when LOH occurring and not occurring.

Using `dwgsim`, I first 
```
dwgsim -d 300 -s 20 -N 10000 -1 150 -2 150 -R 0.02 -S 2 -H -o 1 -e 0.00109 -E 0.00109 s6.n.hla.fasta s6.n
```

```
dwgsim -d 280 -s 20 -N 20000 -1 150 -2 150 -R 0.02 -S 2 -H -o 1 -e 0.002 -E 0.002 s6.t.hla.fasta s6.t
```

Simulate the lost allele (`hla_a_30_01_01`) with low coverage:
```
dwgsim -d 280 -s 20 -N 200 -1 150 -2 150 -R 0.02 -S 2 -H -o 1 -e 0.002 -E 0.002 s6.t.a2.fasta s6.t.a2

```
