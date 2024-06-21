## Simulation dataset

I use [dwgsim](https://github.com/nh13/DWGSIM/tree/main) to simulate a few datasets as guidance set up expectation when LOH occurring and not occurring.

## Case: Subject s6

Subject `s6` loses an allele of HLA-A gene, and an amplification of one allele at HLA-C gene locus. The copy number (CN) is coded as 1/0; 1/1; 3/1.

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

```
hlalohReforged --tbam s6_t.hla.realn.ready.bam --nbam s6_n.hla.realn.ready.bam s6_n.hla.fasta --tstates tstates.tsv --outdir $PWD --subject s6 --min_necnt 2 --corrector global
```

## Case: copy-neutral LOH

## Case: Non-LOH


## Case: 1/1;2/1;0/0

Simulate a case where HLA-C gene locus is entirely deleted
```
hlalohReforged --tbam s7_t.hla.realn.ready.bam --nbam s7_n.hla.realn.ready.bam s7_n.hla.fasta --tstates tstates.tsv --outdir $PWD --subject s7 --min_necnt 2 --corrector global
```

