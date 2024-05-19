require(argparse)

parse_cmd <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--subject",
    metavar = "STR", type = "character", required = TRUE,
    help = "Specify the subject ID"
  )
  parser$add_argument("--tbam",
    metavar = "FILE", type = "character", required = TRUE,
    help = "Specify the tumor bam file"
  )
  parser$add_argument("--nbam",
    metavar = "FILE", type = "character", required = TRUE,
    help = "Specify the normal bam file"
  )
  parser$add_argument("--hlaref",
    metavar = "FILE", type = "character", required = TRUE,
    help = "Specify HLA reference sequence"
  )
  parser$add_argument("--tstates",
    metavar = "FILE", type = "character",
    help = "Specify file includeing tumor purity and ploidy"
  )
  parser$add_argument("--outdir",
    metavar = "DIR", type = "character", required = TRUE,
    help = "Specify the output directory"
  )
  parser$add_argument("--min_cov",
    metavar = "INT", type = "integer", default = 30,
    help = "Specify the minimum coverage at mismatch sites (30)"
  )
  parser$add_argument("--min_necnt",
    metavar = "INT", type = "integer", default = 1,
    help = paste(
      "Specify the minimum number of diff events",
      "allowed for reads mapping to HLA alleles (1)"
    )
  )
  parser$add_argument("--threads",
    metavar = "INT", type = "integer", default = 16,
    help = "Specify the number of threads (16)"
  )
  # parser$add_argument("--example",
  #  action = "store_true",
  #  help = "Specify to run on example data provided by LOHHLA"
  # )

  parser$parse_args()
}
