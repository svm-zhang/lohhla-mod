#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(ggplot2)
  library(splitstackshape)
  library(seqinr)
  library(Biostrings)
  library(Rsamtools)
})

options(width = 600)


extract_tstates <- function(tstate_est_file) {
  print("[INFO] Getting estimates of ploidy and purity")
  tstates_dt <- fread(tstate_est_file, drop = c(1))
  if (!all(c("TumorPloidy", "TumorPurityNGS") %in% names(tstates_dt))) {
    print(paste(
      "[ERROR] Miss either purity, or ploidy, or both ",
      "in the CNV model results",
      sep = ""
    ))
    quit(status = 1)
  }
  list(
    ploidy = tstates_dt[["TumorPloidy"]],
    purity = tstates_dt[["TumorPurityNGS"]]
  )
}

make_bins_old <- function(allele, aln, allele_length, bin_size = 150) {
  start_pos <- aln$start
  end_pos <- aln$end
  bin_breaks <- seq(start_pos, end_pos, by = bin_size)
  # the last bin can be less than 150, so merged with second last bin
  # end_pos + 2 was from original code
  bin_breaks <- c(bin_breaks[-length(bin_breaks)], end_pos + 2)
  istarts <- bin_breaks[-length(bin_breaks)]
  # this makes sure no end position of last bin not repeating as
  # start position in the next bin
  istarts[2:length(istarts)] <- istarts[2:length(istarts)] + 1
  iends <- bin_breaks[2:length(bin_breaks)]
  indices <- seq(1, length(istarts))

  if (start_pos > 2) {
    istarts <- c(1, istarts)
    iends <- c(start_pos - 1, iends)
    indices <- c(0, indices)
  }
  if (allele_length > max(iends)) {
    istarts <- c(istarts, max(iends) + 1)
    iends <- c(iends, allele_length)
    indices <- c(indices, allele_length + 1)
  }

  bin_dt <- data.table(
    seqnames = allele, start = istarts, end = iends, bin = indices
  )
  bin_dt[, end := ifelse(end > allele_length, allele_length, end)]
  print(bin_dt)
  stop()

  bin_dt
}

initialize_libs <- function() {
  pkg_name <- "lohhlamod"
  libpaths <- file.path(.libPaths(), pkg_name)
  lib_pattern <- "(bamer|cli|loh|pairwise_aln|pathio).R"
  rscripts <- list.files(libpaths, lib_pattern, full.names = TRUE, recursive = TRUE)
  if (length(rscripts) == 0) {
    print("[ERROR] Found no related R script libraries to run lohhlamod")
    quit(status = 1)
  }
  invisible(sapply(rscripts, source))
}

main <- function() {
  initialize_libs()

  gamma <- 1

  args <- parse_cmd()

  parse_file_path(file = args$tbam)
  parse_file_path(file = args$nbam)
  parse_dir_path(dir = args$outdir, create = TRUE)

  tid <- extract_rgsm_from_bam_header(bam = args$tbam)
  nid <- extract_rgsm_from_bam_header(bam = args$nbam)

  purity <- 0.5
  ploidy <- 2
  if (!is.null(args$tstates)) {
    tstates <- extract_tstates(tstate_est_file = args$tstates)
    if (!is.na(tstates$ploidy) && !is.na(tstates$purity)) {
      ploidy <- tstates$ploidy
      purity <- tstates$purity
    }
  }
  print(paste("[INFO] Purity = ", purity, " Ploidy = ", ploidy, sep = ""))

  alleles_n <- extract_seqinfo_from_bam(bam = args$nbam)
  alleles_t <- extract_seqinfo_from_bam(bam = args$tbam)
  # extract function returns named integer vector
  # here we need names of that vector to get allele name
  alleles_n <- sort(names(alleles_n))
  alleles_t <- sort(names(alleles_t))
  if (!all.equal(alleles_n, alleles_t)) {
    print("[ERROR] Different set of alleles detected in normal and tumor BAMs")
    print("[ERROR] Alleles in normal BAM: ", paste(alleles_n, collapse = "|"))
    print("[ERROR] Alleles in tumor BAM: ", paste(alleles_t, collapse = "|"))
    quit(status = 1)
  }

  alleles_dt <- data.table(Alleles = alleles_n)
  alleles_dt[, HLAGene := tstrsplit(Alleles, "_", keep = 2)]
  alleles_dt[, HLAGene := paste("hla_", HLAGene, sep = "")]
  alleles_dt[, Alleles := paste(Alleles, collapse = ","), by = HLAGene]
  alleles_dt[, c("A1", "A2") := tstrsplit(Alleles, ","), by = HLAGene]
  alleles_dt[, Alleles := NULL]
  alleles_dt <- unique(alleles_dt, by = "HLAGene")

  # filter input bams by ecnt
  filt_nbam <- file.path(args$outdir, paste(nid, ".filt.bam", sep = ""))
  filt_tbam <- file.path(args$outdir, paste(tid, ".filt.bam", sep = ""))
  if (!file.exists(filt_nbam) || !file.exists(filt_tbam)) {
    filter_bam_by_ecnt(
      bam = args$nbam, obam = filt_nbam, min_necnt = args$min_necnt
    )
    filter_bam_by_ecnt(
      bam = args$tbam, obam = filt_tbam, min_necnt = args$min_necnt
    )
  }

  n_seq_depth <- estimate_dp(filt_nbam, alleles = alleles_n)
  t_seq_depth <- estimate_dp(filt_tbam, alleles = alleles_n)
  corrector <- n_seq_depth / t_seq_depth
  print(paste("[INFO] Global depth corrector: ", corrector, sep=""))

  loh_res_dt <- alleles_dt[, call_hla_loh(
    .SD,
    tbam = filt_tbam, nbam = filt_nbam, hlaref = args$hlaref,
    outdir = args$outdir, purity = purity, ploidy = ploidy,
    min_dp = args$min_cov, min_necnt = args$min_nm,
    corrector = corrector, gamma = gamma
  ), by = "HLAGene"]
  out_res <- file.path(args$outdir, paste(args$subject, ".loh.res.tsv", sep = ""))
  fwrite(loh_res_dt, out_res, sep = "\t", row.names = FALSE, quote = FALSE)
}

main()
