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
  parser$add_argument("--example",
    action = "store_true",
    help = "Specify to run on example data provided by LOHHLA"
  )

  parser$parse_args()
}

t_test_with_na <- function(x, alternative = "two.sided", mu = 0) {
  if (length(x[!is.na(x)]) <= 1) {
    stat <- NA
    ci <- c(NA, NA)
    out <- list(stat, ci)
    names(out) <- c("stat", "conf.int")
    out
  } else {
    t.test(x, alternative = alternative, mu = mu)
  }
}

test_cn_loss <- function(x) {
  if (anyNA(x)) {
    return(NA)
  }
  x <- matrix(round(x), nrow = 2)
  test <- fisher.test(x)

  test$p.value
}

parse_file_path <- function(file) {
  if (!file.exists(file)) {
    print(paste("[ERROR] Cannot find the file provided: ", file, sep = ""))
    quit(status = 1)
  }
  file <- normalizePath(file, mustWork = TRUE)
}

parse_dir_path <- function(dir, create) {
  if (!dir.exists(dir)) {
    if (create != TRUE) {
      print(paste("[ERROR] Cannot find the file provided: ", dir, sep = ""))
    } else {
      dir.create(dir, recursive = TRUE)
      dir <- normalizePath(dir, mustWork = TRUE)
    }
  }
}

extract_bam_header <- function(bam) {
  bf <- BamFile(bam)
  scanBamHeader(bf)
}

extract_seqinfo_from_bam <- function(bam) {
  header <- extract_bam_header(bam = bam)
  # a named integer vector
  targets <- header$targets
  targets
}

extract_rg_from_bam_header <- function(bam) {
  header <- extract_bam_header(bam = bam)
  rg <- header$text$`@RG`
  if (is.null(rg) || is.na(rg)) {
    print(paste("[ERROR] Cannot find the read group (@RG) in the BAM file: ",
      bam,
      sep = ""
    ))
    quit(status = 1)
  }
  rg
}

extract_rgsm_from_bam_header <- function(bam) {
  rg <- extract_rg_from_bam_header(bam = bam)
  sm <- rg[which(grepl("^SM", rg))]
  if (length(sm) == 0) {
    print(paste("[ERROR] Cannot find the SM in the @RG: ",
      paste(rg, collapse = "\t"),
      sep = ""
    ))
    quit(status = 1)
  }
  if (length(sm) > 1) {
    print("[INFO] Expect to extract only one read group in BAM")
    print("[INFO] Only the first one will be extracted")
    sm <- sm[1]
  }
  sm <- gsub("^SM:", "", sm)
  sm
}

paste_vector <- function(v, sep = "") {
  vt <- v[1]
  if (length(v) > 1) {
    for (g in 2:length(v)) {
      vt <- paste(vt, v[g], sep = sep)
    }
  }
  vt <- paste(vt, " EnD", sep = "")
  out_v <- sub(" EnD", "", vt)
  out_v <- sub("NA , ", "", out_v)
  out_v <- sub(" , NA", "", out_v)
  out_v <- sub(" , NA , ", " , ", out_v)
  out_v
}

get_mm_bw_alleles <- function(alignment, chunksize = 60, returnlist = FALSE) {
  a1_aln <- pattern(alignment) # Get the alignment for the first sequence
  a2_aln <- subject(alignment) # Get the alignment for the second sequence

  a1_aln_start <- start(pattern(alignment))
  a1_aln_end <- end(pattern(alignment))
  a2_aln_start <- start(subject(alignment))
  a2_aln_end <- end(subject(alignment))

  k <- 1 + a1_aln_start - 1
  # a1_aln_seq is a character vector showing the alignment from a1 pov
  a1_aln_seq <- unlist(strsplit(as.character(a1_aln), split = ""))
  seq1_positions <- c()
  for (char in a1_aln_seq) {
    if (char %in% c("C", "G", "A", "T")) {
      seq1_positions <- c(seq1_positions, k)
      k <- k + 1
      next
    }

    if (char %in% c("-")) {
      seq1_positions <- c(seq1_positions, k)
      next
    }
  }

  k <- 1 + a2_aln_start - 1
  a2_aln_seq <- unlist(strsplit(as.character(a2_aln), split = ""))
  seq2_positions <- c()
  for (char in a2_aln_seq) {
    if (char %in% c("C", "G", "A", "T")) {
      seq2_positions <- c(seq2_positions, k)
      k <- k + 1
      next
    }

    if (char %in% c("-")) {
      seq2_positions <- c(seq2_positions, k)
      next
    }
  }

  seq1_diff <- seq1_positions[a1_aln_seq != a2_aln_seq]
  seq2_diff <- seq2_positions[a1_aln_seq != a2_aln_seq]

  type1_diff <- rep(1, length(seq1_diff))
  type1_diff[which(a1_aln_seq[a1_aln_seq != a2_aln_seq] %in% "-")] <- 2

  type2_diff <- rep(1, length(seq2_diff))
  type2_diff[which(a2_aln_seq[a2_aln_seq != a1_aln_seq] %in% "-")] <- 2

  aln_dt <- data.table(
    seq1_diff = seq1_diff,
    seq2_diff = seq2_diff,
    type1_diff = type1_diff,
    type2_diff = type2_diff
  )
  aln_dt <- aln_dt[type1_diff != 2 & type2_diff != 2]

  out <- list()
  out$diffSeq1 <- aln_dt$seq1_diff
  out$diffSeq2 <- aln_dt$seq2_diff
  out$a1 <- list(start = a1_aln_start, end = a1_aln_end)
  out$a2 <- list(start = a2_aln_start, end = a2_aln_end)

  out
}

get_mismatches_bw_alleles <- function(a1_seq, a2_seq) {
  a1_seq <- paste_vector(toupper(a1_seq), sep = "")
  a2_seq <- paste_vector(toupper(a2_seq), sep = "")
  sigma <- nucleotideSubstitutionMatrix(
    match = 2, mismatch = -1, baseOnly = TRUE
  )
  pair_aln <- pairwiseAlignment(
    a1_seq, a2_seq,
    substitutionMatrix = sigma, gapOpening = -2, gapExtension = -4,
    scoreOnly = FALSE, type = "local"
  )
  mm <- get_mm_bw_alleles(pair_aln, returnlist = TRUE)

  mm
}

get_indel_length <- function(cigar) {
  tmp <- unlist(strsplit(gsub("([0-9]+)", "~\\1~", cigar), "~"))
  ins <- grep(pattern = "I", x = tmp)
  del <- grep(pattern = "D", x = tmp)
  total <- sum(as.numeric(tmp[(ins - 1)])) + sum(as.numeric(tmp[del - 1]))
  total
}

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

count_n_reads_from_bam <- function(
    bam, count_dup = FALSE, only_proper = TRUE, only_primary = TRUE,
    tagfilter = list()) {
  count_n_reads_df <- countBam(
    bam,
    param = ScanBamParam(
      flag = scanBamFlag(
        isDuplicate = count_dup,
        isProperPair = only_proper,
        isSecondaryAlignment = !only_primary
      ),
      tagFilter = tagfilter
    )
  )
  count_n_reads_df$records
}

filter_bam_by_ecnt <- function(bam, obam, min_necnt = 1) {
  bamf <- BamFile(file = bam)

  scan_param <- ScanBamParam(
    flag = scanBamFlag(),
    what = c("qname", "flag", "cigar"),
    tag = "NM"
  )
  aln <- scanBam(bamf, param = scan_param)
  aln_dt <- data.table(
    qname = aln[[1]]$qname,
    cigar = aln[[1]]$cigar,
    flag = aln[[1]]$flag,
    nm = unlist(aln[[1]]$tag)
  )
  aln_dt[, read_idx := ifelse(
    bamFlagAsBitMatrix(as.integer(flag))[7] == 1, 1, 2
  ),
  by = seq_len(nrow(aln_dt))
  ]
  if (nrow(aln_dt) == 0) {
    print(paste(
      "[ERROR] Found no alignments in the give BAM: ", bam,
      sep = ""
    ))
    quit(status = 1)
  }
  cigar <- NULL # this is to avoid "no visible binding for cigar"
  n_ins <- n_del <- n_mm <- n_ecnt <- NULL
  nread_per_frag <- NULL
  aln_dt[, ":="(
    n_ins = length(
      grep(pattern = "I", unlist(strsplit(cigar, "")))
    ),
    n_del = length(
      grep(pattern = "D", unlist(strsplit(cigar, "")))
    )
  ), by = seq_len(nrow(aln_dt))]
  aln_dt[, "n_mm" := nm - apply(.SD, 1, get_indel_length), .SDcols = "cigar"]
  aln_dt[, "n_ecnt" := n_mm + n_ins + n_del]
  aln_dt <- aln_dt[n_ecnt <= min_necnt]
  aln_dt[, "nread_per_frag" := .N, by = "qname"]
  aln_dt <- aln_dt[nread_per_frag == 2]
  filter <- S4Vectors::FilterRules(
    list(function(x) x$qname %in% aln_dt$qname)
  )
  obam <- filterBam(bamf, obam, filter = filter, param = scan_param)
}

init_loh_report <- function(a1, a2) {
  list(
    "HLA_A1" = a1,
    "HLA_A2" = a2,
    "HLA_A1_CN" = NaN,
    "HLA_A1_CN_Lower" = NaN,
    "HLA_A1_CN_Upper" = NaN,
    "HLA_A2_CN" = NaN,
    "HLA_A2_CN_Lower" = NaN,
    "HLA_A2_CN_Upper" = NaN,
    "HLA_A1_Median_LogR" = NaN,
    "HLA_A2_Median_LogR" = NaN,
    "HLA_A1_MM_Median_LogR" = NaN,
    "HLA_A2_MM_Median_LogR" = NaN,
    "MM_LogR_Paired_Pvalue" = NaN,
    "Median_BAF" = NaN,
    "Num_MM" = as.integer(0),
    "Num_Bins" = as.integer(0),
    "Num_CN_Loss_Supporting_Bins" = as.integer(0)
  )
}

get_allele_coverage <- function(allele, bam, min_dp = 0) {
  print(paste(
    "[INFO] Get coverage for ", allele, " from ", bam,
    sep = ""
  ))
  seqinfo <- extract_seqinfo_from_bam(bam = bam)

  allele_seq_ln <- seqinfo[which(names(seqinfo) == allele)]
  allele_to_scan <- GenomicRanges::GRanges(
    seqnames = allele,
    ranges = IRanges::IRanges(start = 1, end = allele_seq_ln)
  )
  scanflag <- scanBamFlag(isProperPair = TRUE)
  scan_param <- ScanBamParam(
    flag = scanflag,
    what = c("qname", "flag"),
    which = allele_to_scan,
  )
  pileup_param <- PileupParam(
    min_mapq = 20, distinguish_strands = FALSE
  )
  p_dt <- setDT(
    pileup(
      file = bam, scanBamParam = scan_param, pileupParam = pileup_param
    )
  )
  # pileup return seqnames as factor
  p_dt <- p_dt[count > min_dp]
  p_dt[, seqnames := as.character(seqnames)]
  p_dt[, which_label := NULL]
  p_dt
}

prep_allelic_cov <- function(t_dt, n_dt, bin_dt, multfactor) {
  cov_dt <- combine_tn_cov(t_dt = t_dt, n_dt = n_dt)
  cov_dt <- bin_allele_cov(cov_dt = cov_dt, bin_dt = bin_dt)
  cov_dt <- estimate_logr(cov_dt = cov_dt, multfactor = multfactor)
  cov_dt
}

combine_tn_cov <- function(t_dt, n_dt) {
  a_t <- unique(t_dt$seqnames)
  a_n <- unique(n_dt$seqnames)
  if (length(a_t) != 1 || length(a_n) != 1) {
    print("[ERROR] Only one seqname is expected in tumor and normal tables")
    quit(status = 1)
  }
  if (a_t != a_n) {
    print("[ERROR] To be combined tumor and normal tables have diff seqnames")
    quit(status = 1)
  }
  setkey(t_dt, seqnames, pos, nucleotide)
  setkey(n_dt, seqnames, pos, nucleotide)
  allele_cov_dt <- t_dt[n_dt]
  allele_cov_dt[, ":="(t_dp = count, n_dp = i.count)]
  allele_cov_dt[, ":="(count = NULL, i.count = NULL)]

  allele_cov_dt
}

make_bins <- function(allele, aln, allele_length, bin_size = 150) {
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

  bin_dt
}

bin_allele_cov <- function(cov_dt, bin_dt) {
  # make granges for allele coverage dt and bins
  tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(
    cov_dt,
    start.field = "pos", end.field = "pos"
  )
  bin_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bin_dt,
    keep.extra.columns = TRUE
  )
  ovl <- GenomicRanges::findOverlaps(tmp_gr, bin_gr)
  cov_dt[S4Vectors::queryHits(ovl), "bin"] <- bin_dt[
    S4Vectors::subjectHits(ovl), "bin"
  ]
  cov_dt[, ":="(
    bin_t_dp = as.numeric(median(t_dp, na.rm = TRUE)),
    bin_n_dp = as.numeric(median(n_dp, na.rm = TRUE))
  ), by = "bin"]
  cov_dt
}

estimate_logr <- function(cov_dt, multfactor) {
  cov_dt[, logR := log2(t_dp / n_dp * multfactor)]
  # FIXME: i can also use bin_t_dp and bin_n_dp calculate
  # two should be similar
  cov_dt[, bin_logR := median(logR, na.rm = TRUE), by = "bin"]
  cov_dt
}

estimate_binned_logr <- function(a1_dt, a2_dt, multfactor) {
  bin_dt <- merge(
    unique(a1_dt, by = "a1_bin"),
    unique(a2_dt, by = "a2_bin"),
    by.x = c("a1_bin"), by.y = c("a2_bin"),
    all.x = TRUE, all.y = TRUE
  )
  binned_names <- names(bin_dt)[which(grepl("bin", names(bin_dt)))]
  bin_dt <- bin_dt[, ..binned_names]
  bin_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
  bin_dt[, logR_combined_bin := log2(
    (a1_bin_t_dp + a2_bin_t_dp) /
      (a1_bin_n_dp + a2_bin_n_dp) * multfactor
  )]
  bin_dt
}

prep_mm_cov <- function(mm, a1_dt, a2_dt) {
  mm_est_dt <- data.table(a1_pos = mm$diffSeq1, a2_pos = mm$diffSeq2)
  mm_est_dt <- mm_est_dt[
    a1_pos %in% a1_dt$a1_pos & a2_pos %in% a2_dt$a2_pos
  ]
  if (nrow(mm_est_dt) == 0) {
    return(NULL)
  }
  simply_cols <- names(a1_dt)[which(!grepl("bin", names(a1_dt)))]
  simply_cols <- c(simply_cols, "a1_bin")
  mm_est_dt <- merge(
    mm_est_dt, a1_dt[, ..simply_cols],
    by = "a1_pos", all.x = TRUE
  )
  mm_est_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
  if (nrow(mm_est_dt) == 0) {
    return(NULL)
  }
  simply_cols <- names(a2_dt)[which(!grepl("bin", names(a2_dt)))]
  mm_est_dt <- merge(
    mm_est_dt, a2_dt[, ..simply_cols],
    by = "a2_pos", all.x = TRUE
  )
  if (nrow(mm_est_dt) == 0) {
    return(NULL)
  }

  mm_est_dt
}

estimate_baf <- function(mm_dt, bin_dt) {
  mm_dt[, baf := a1_t_dp / (a1_t_dp + a2_t_dp)] # nolint
  if (!"capture_bias_bin" %in% names(bin_dt)) {
    print("[WARN] Found no column named capture_bias_bin in bin_dt")
    print("[WARN] No capture bias will be corrected for BAF")
    mm_dt[, baf_correct := baf]
    return(mm_dt)
  }
  req_cols <- c("bin", "logR_combined_bin")
  miss_cols <- req_cols[which(!req_cols %in% names(bin_dt))]
  if (length(miss_cols) > 0) {
    print(paste(
      "[ERROR] Miss ",
      paste(miss_cols, collapse = ","),
      " columns in bin_dt",
      sep = ""
    ))
    print("[ERROR] Cannot continue estimate BAF")
    quit(status = 1)
  }
  mm_dt <- merge(
    mm_dt,
    bin_dt[, c("bin", "logR_combined_bin", "capture_bias_bin")],
    by = "bin",
    all.x = TRUE
  )
  mm_dt[, baf_correct := baf / capture_bias_bin]
  mm_dt
}

estimate_cn <- function(mm_dt) {
  mm_dt[, "a1_cn" :=
    (purity - 1 + baf_correct * 2^(logR_combined_bin / gamma) * # nolint
      ((1 - purity) * 2 + purity * ploidy)) / purity, ]
  mm_dt[, "a2_cn" :=
    (purity - 1 - (baf_correct - 1) * 2^(logR_combined_bin / gamma) * # nolint
      ((1 - purity) * 2 + purity * ploidy)) / purity, ]
  mm_dt
}

estimate_cn_conf <- function(cn_dt, which) {
  col <- NULL
  pattern <- paste(which, "bin_cn", sep = "_")
  col <- names(cn_dt)[which(grepl(pattern, names(cn_dt)))]
  if (is.null(col) || length(col) == 0) {
    print(paste(
      "[ERROR] Failed to find cn column in cn_dt for ",
      which, " allele",
      sep = ""
    ))
    quit(status = 1)
  }
  cn_test <- t_test_with_na(cn_dt[[col]])
  cn_est_conf <- cn_test$conf.int
  cn_est_lower <- cn_est_conf[1]
  cn_est_upper <- cn_est_conf[2]
  list(est_lower = cn_est_lower, est_upper = cn_est_upper)
}

dump_intermedia_tables <- function(out, mm, cov_dt, bin_dt = NULL, mm_dt = NULL) {
  to_save <- list(
    mm = mm,
    cov_dt = cov_dt,
    bin_dt = bin_dt,
    mm_dt = mm_dt
  )
  saveRDS(to_save, out)
}

call_hla_loh <- function(
    dt, tbam, nbam, hlaref, outdir,
    purity, ploidy, multfactor, min_dp, min_necnt,
    tid = "example_tumor", nid = "example_normal", gamma = 1) {
  a1 <- dt$A1
  a2 <- dt$A2
  alleles_str <- paste(c(a1, a2), collapse = " and ")
  print(paste("[INFO] Analyze LOH for ", alleles_str, sep = ""))

  report <- init_loh_report(a1, a2)

  hla_seq <- read.fasta(hlaref)
  a1_seq <- hla_seq[[a1]]
  a2_seq <- hla_seq[[a2]]
  print(paste(
    "[INFO] Align sequences between ", alleles_str,
    sep = ""
  ))
  mm <- get_mismatches_bw_alleles(a1_seq, a2_seq)
  report$Num_MM <- length(mm$diffSeq1)

  if (length(mm$diffSeq1) == 0) {
    print(paste(
      "[INFO] there is no difference between the two alleles",
      a1, a2,
      sep = " "
    ))
    print("[INFO] no call will be made. Move to the next HLA gene")
    return(report)
  }

  if (length(mm$diffSeq1) < 5) {
    print("[INFO] HLA alleles are similar (less than 5 mismatch positions)")
    print("[INFO] no call will be made. Move to the next HLA gene")
    return(report)
  }

  print("[INFO] Make bins accounting for alignment start position")
  a1_bin_dt <- make_bins(
    allele = a1,
    aln = mm$a1,
    allele_length = length(a1_seq)
  )
  a2_bin_dt <- make_bins(
    allele = a2,
    aln = mm$a2,
    allele_length = length(a2_seq)
  )

  t_a1_cov <- get_allele_coverage(allele = a1, bam = tbam)
  t_a2_cov <- get_allele_coverage(allele = a2, bam = tbam)
  n_a1_cov <- get_allele_coverage(allele = a1, bam = nbam, min_dp = min_dp)
  n_a2_cov <- get_allele_coverage(allele = a2, bam = nbam, min_dp = min_dp)
  t_a1_cov <- t_a1_cov[pos %in% n_a1_cov$pos]
  t_a2_cov <- t_a2_cov[pos %in% n_a2_cov$pos]
  if (nrow(n_a1_cov) == 0 || nrow(n_a2_cov) == 0) {
    print("[INFO] Found no coverage for either allele in normal")
    print("[INFO] Move to next HLA gene")
    return(report)
  }
  print(paste(
    "[INFO] Prepare coverage table for allele ", a1,
    sep = ""
  ))
  a1_cov_dt <- prep_allelic_cov(
    t_dt = t_a1_cov, n_dt = n_a1_cov,
    bin_dt = a1_bin_dt, multfactor = multfactor
  )
  names(a1_cov_dt) <- paste("a1_", names(a1_cov_dt), sep = "")
  print(paste(
    "[INFO] Prepare coverage table for allele ", a2,
    sep = ""
  ))
  a2_cov_dt <- prep_allelic_cov(
    t_dt = t_a2_cov, n_dt = n_a2_cov,
    bin_dt = a2_bin_dt, multfactor = multfactor
  )
  names(a2_cov_dt) <- paste("a2_", names(a2_cov_dt), sep = "")
  report$HLA_A1_Median_LogR <- median(a1_cov_dt$a1_logR, na.rm = TRUE)
  report$HLA_A2_Median_LogR <- median(a2_cov_dt$a2_logR, na.rm = TRUE)
  report$HLA_A1_MM_Median_LogR <- median(
    a1_cov_dt[a1_pos %in% mm$diffSeq1]$a1_logR,
    na.rm = TRUE
  )
  report$HLA_A2_MM_Median_LogR <- median(
    a2_cov_dt[a2_pos %in% mm$diffSeq2]$a2_logR,
    na.rm = TRUE
  )

  bin_dt <- estimate_binned_logr(
    a1_dt = a1_cov_dt, a2_dt = a2_cov_dt, multfactor = multfactor
  )
  bin_dt[, cn_loss_test_bin := apply(
    .SD, 1, test_cn_loss
  ),
  .SDcols = c("a1_bin_t_dp", "a1_bin_n_dp", "a2_bin_t_dp", "a2_bin_n_dp")
  ]
  bin_dt[, capture_bias_bin := a1_bin_n_dp / a2_bin_n_dp]
  setkey(bin_dt, bin)
  report$Num_Bins <- nrow(bin_dt)
  report$Num_CN_Loss_Supporting_Bins <- nrow(
    bin_dt[cn_loss_test_bin <= 0.01]
  )

  mm_dt <- prep_mm_cov(mm = mm, a1_dt = a1_cov_dt, a2_dt = a2_cov_dt)
  if (is.null(mm_dt)) {
    print("[INFO] No coverage info found for both alleles at mm sites")
    print("[INFO] Cannot estimate copy numbers. Moving to next HLA gene")
    return(report)
  }
  mm_dt <- estimate_baf(mm_dt = mm_dt, bin_dt = bin_dt)
  paired_t_test <- t.test(mm_dt$a1_logR, mm_dt$a2_logR, paired = TRUE)
  report$MM_LogR_Paired_Pvalue <- paired_t_test$p.value
  report$Median_BAF <- median(mm_dt$baf_correct, na.rm = TRUE)

  print("[INFO] Estimate copy number at mismatch sites")
  mm_dt <- estimate_cn(mm_dt = mm_dt)
  mm_dt[, ":="(
    a1_bin_cn = median(a1_cn, na.rm = TRUE),
    a2_bin_cn = median(a2_cn, na.rm = TRUE)
  ),
  by = "bin"
  ]
  cn_dt <- unique(mm_dt, by = "bin")
  setkey(cn_dt, bin)
  report$HLA_A1_CN <- round(median(cn_dt$a1_bin_cn, na.rm = TRUE), 4)
  report$HLA_A2_CN <- round(median(cn_dt$a2_bin_cn, na.rm = TRUE), 4)

  cn_est_conf <- estimate_cn_conf(cn_dt = cn_dt, which = "a1")
  report$HLA_A1_CN_Lower <- round(cn_est_conf$est_lower, 4)
  report$HLA_A1_CN_Upper <- round(cn_est_conf$est_upper, 4)
  cn_est_conf <- estimate_cn_conf(cn_dt = cn_dt, which = "a2")
  report$HLA_A2_CN_Lower <- round(cn_est_conf$est_lower, 4)
  report$HLA_A2_CN_Upper <- round(cn_est_conf$est_upper, 4)
  bin_dt <- cn_dt[, c("bin", "a1_bin_cn", "a2_bin_cn")][bin_dt]
  rm(cn_dt)

  hla_gene <- gsub("(_|\\*)*(([0-9])+(_|:)*)+", "", a1)
  hla_gene <- tolower(hla_gene)
  out_rds <- file.path(outdir, paste(hla_gene, ".rds", sep = ""))
  names(a1_cov_dt) <- gsub("a1_", "", names(a1_cov_dt))
  names(a2_cov_dt) <- gsub("a2_", "", names(a2_cov_dt))
  dump_intermedia_tables(
    out = out_rds,
    mm = mm,
    cov_dt = rbind(a1_cov_dt, a2_cov_dt),
    bin_dt = bin_dt,
    mm_dt = mm_dt
  )


  report
}

gamma <- 1
bin_size <- 150

args <- parse_cmd()

parse_file_path(file = args$tbam)
parse_file_path(file = args$nbam)
parse_dir_path(dir = args$outdir, create = TRUE)

tid <- NULL
nid <- NULL
if (args$example) {
  # when using the test data provided by LOHHLA.R set these two lines below
  tid <- "example_tumor"
  nid <- "example_normal"
} else {
  tid <- extract_rgsm_from_bam_header(bam = args$tbam)
  nid <- extract_rgsm_from_bam_header(bam = args$nbam)
}

if (args$example) {
  # when using the test data provided by LOHHLA.R set these two lines below
  tid <- "example_tumor"
  nid <- "example_normal"
}

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

print("[INFO] Counting sequencing depth from realigned normal BAM")
# FIXME: tagfilter should be generated with args$min_necnt, rather than
# hard-coded
n_seq_depth <- count_n_reads_from_bam(
  bam = args$nbam, tagfilter = list(NM = c(0, 1))
)
print("[INFO] Counting sequencing depth from realigned tumor BAM")
t_seq_depth <- count_n_reads_from_bam(
  bam = args$tbam, tagfilter = list(NM = c(0, 1))
)
if (t_seq_depth <= 0) {
  print("[ERROR] Found no alignments in the provided tumor BAM")
  quit(status = 1)
}
print(n_seq_depth)
print(t_seq_depth)
multfactor <- n_seq_depth / t_seq_depth
print(paste("[INFO] Normal sequencing depth = ", n_seq_depth, sep = ""))
print(paste("[INFO] Tumor sequencing depth = ", t_seq_depth, sep = ""))
print(paste(
  "[INFO] Sequencing depth correcting factor = ", multfactor,
  sep = ""
))

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

loh_res_dt <- alleles_dt[, call_hla_loh(
  .SD,
  tbam = filt_tbam, nbam = filt_nbam, hlaref = args$hlaref,
  outdir = args$outdir, purity = purity, ploidy = ploidy,
  multfactor = multfactor, min_dp = args$min_cov, min_necnt = args$min_nm,
  tid = tid, nid = nid
), by = "HLAGene"]
out_res <- file.path(args$outdir, paste(tid, ".loh.res.tsv", sep = ""))
fwrite(loh_res_dt, out_res, sep = "\t", row.names = FALSE, quote = FALSE)
