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

# parse_cmd <- function() {
#  parser <- ArgumentParser()
#  parser$add_argument("--subject",
#    metavar = "STR", type = "character", required = TRUE,
#    help = "Specify the subject ID"
#  )
#  parser$add_argument("--tbam",
#    metavar = "FILE", type = "character", required = TRUE,
#    help = "Specify the tumor bam file"
#  )
#  parser$add_argument("--nbam",
#    metavar = "FILE", type = "character", required = TRUE,
#    help = "Specify the normal bam file"
#  )
#  parser$add_argument("--hlaref",
#    metavar = "FILE", type = "character", required = TRUE,
#    help = "Specify HLA reference sequence"
#  )
#  parser$add_argument("--tstates",
#    metavar = "FILE", type = "character",
#    help = "Specify file includeing tumor purity and ploidy"
#  )
#  parser$add_argument("--outdir",
#    metavar = "DIR", type = "character", required = TRUE,
#    help = "Specify the output directory"
#  )
#  parser$add_argument("--min_cov",
#    metavar = "INT", type = "integer", default = 30,
#    help = "Specify the minimum coverage at mismatch sites (30)"
#  )
#  parser$add_argument("--min_necnt",
#    metavar = "INT", type = "integer", default = 1,
#    help = paste(
#      "Specify the minimum number of diff events",
#      "allowed for reads mapping to HLA alleles (1)"
#    )
#  )
#  parser$add_argument("--threads",
#    metavar = "INT", type = "integer", default = 16,
#    help = "Specify the number of threads (16)"
#  )
#  parser$add_argument("--example",
#    action = "store_true",
#    help = "Specify to run on example data provided by LOHHLA"
#  )
#
#  parser$parse_args()
# }

# t_test_with_na <- function(x, alternative = "two.sided", mu = 0) {
#  if (length(x[!is.na(x)]) <= 1) {
#    stat <- NA
#    ci <- c(NA, NA)
#    out <- list(stat, ci)
#    names(out) <- c("stat", "conf.int")
#    out
#  } else {
#    t.test(x, alternative = alternative, mu = mu)
#  }
# }

# test_cn_loss <- function(x) {
#  if (anyNA(x)) {
#    return(NA)
#  }
#  x <- matrix(round(x), nrow = 2)
#  test <- fisher.test(x)
#
#  test$p.value
# }

# parse_file_path <- function(file) {
#  if (!file.exists(file)) {
#    print(paste("[ERROR] Cannot find the file provided: ", file, sep = ""))
#    quit(status = 1)
#  }
#  file <- normalizePath(file, mustWork = TRUE)
# }
#
# parse_dir_path <- function(dir, create) {
#  if (!dir.exists(dir)) {
#    if (create != TRUE) {
#      print(paste("[ERROR] Cannot find the file provided: ", dir, sep = ""))
#    } else {
#      dir.create(dir, recursive = TRUE)
#      dir <- normalizePath(dir, mustWork = TRUE)
#    }
#  }
# }

# https://rdrr.io/bioc/Biostrings/src/R/PairwiseAlignments-io.R
# using .pre2postaligned function in writePariwiseAlignments function
# extract_aln_in_pos <- function(axset) {
#  pos <- seq(start(axset@range), end(axset@range))
#  data.table(
#    pos = pos,
#    pa_pos = Biostrings:::.pre2postaligned(pos, axset)
#  )
# }

# get_mismatches_bw_alleles <- function(a1, a2, hlaref) {
#  hla_seq <- read.fasta(hlaref)
#  a1_seq <- hla_seq[[a1]]
#  a2_seq <- hla_seq[[a2]]
#  sigma <- nucleotideSubstitutionMatrix(
#    match = 2, mismatch = -1, baseOnly = TRUE
#  )
#  pair_aln <- pairwiseAlignment(
#    DNAString(paste(toupper(a1_seq), collapse = "")),
#    DNAString(paste(toupper(a2_seq), collapse = "")),
#    substitutionMatrix = sigma, gapOpening = -2, gapExtension = -4,
#    scoreOnly = FALSE, type = "local"
#  )
#  a1_aln <- pattern(pair_aln) # Get the pair_aln for the first sequence
#  a2_aln <- subject(pair_aln) # Get the pair_aln for the second sequence
#
#  a1_aln_start <- start(pattern(pair_aln))
#  a1_aln_end <- end(pattern(pair_aln))
#  a2_aln_start <- start(subject(pair_aln))
#  a2_aln_end <- end(subject(pair_aln))
#
#  p_aln <- extract_aln_in_pos(axset = a1_aln)
#  s_aln <- extract_aln_in_pos(axset = a2_aln)
#  p_aln <- merge(p_aln, s_aln, by = "pa_pos", all = TRUE)
#  p_aln[, ":="(a1_seqnames = a1, a2_seqnames = a2)]
#  p_aln[, a1_ref := unlist(strsplit(as.character(a1_aln), split = ""))]
#  p_aln[, a2_ref := unlist(strsplit(as.character(a2_aln), split = ""))]
#  p_aln <- p_aln[a1_ref != "-" & a2_ref != "-"]
#  setnames(p_aln, c("pos.x", "pos.y"), c("a1_pos", "a2_pos"))
#  diffSeq1 <- p_aln[a1_ref != a2_ref]$a1_pos
#  diffSeq2 <- p_aln[a1_ref != a2_ref]$a2_pos
#
#  if (length(diffSeq1) <= 5) {
#    print("[INFO] Found insufficient num mismatches between the two alleles")
#    return(NULL)
#  }
#  list(
#    aln_dt = p_aln,
#    diffSeq1 = diffSeq1,
#    diffSeq2 = diffSeq2,
#    a1 = list(start = a1_aln_start, end = a1_aln_end, len = length(a1_seq)),
#    a2 = list(start = a2_aln_start, end = a2_aln_end, len = length(a2_seq))
#  )
# }

# get_indel_length <- function(cigar) {
#  tmp <- unlist(strsplit(gsub("([0-9]+)", "~\\1~", cigar), "~"))
#  ins <- grep(pattern = "I", x = tmp)
#  del <- grep(pattern = "D", x = tmp)
#  total <- sum(as.numeric(tmp[(ins - 1)])) + sum(as.numeric(tmp[del - 1]))
#  total
# }

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

# count_n_reads_from_bam <- function(bam, which, tagfilter = list()) {
#  seqinfo <- extract_seqinfo_from_bam(bam = bam)
#  allele_seq_ln <- seqinfo[which(names(seqinfo) %in% which)]
#  dt <- data.table(seqnames = names(allele_seq_ln), end = allele_seq_ln)
#  dt[, start := 1]
#  allele_to_scan <- GenomicRanges::makeGRangesFromDataFrame(dt)
#  count_n_reads_df <- countBam(
#    bam,
#    param = ScanBamParam(
#      which = allele_to_scan,
#      flag = scanBamFlag(
#        isDuplicate = FALSE,
#        isProperPair = TRUE,
#        isNotPassingQualityControls = FALSE,
#        isSupplementaryAlignment = FALSE,
#        isSecondaryAlignment = FALSE
#      ),
#      tagFilter = tagfilter
#    )
#  )
#  sum(count_n_reads_df$records)
# }

# filter_bam_by_ecnt <- function(bam, obam, min_necnt = 1) {
#  bamf <- BamFile(file = bam)
#
#  # including secondary alignments, otherwise, one random
#  # allele out of 2 has more coverage than expected, while
#  # the other has lower-than-expected coverage
#  scanflag <- scanBamFlag(
#    isProperPair = TRUE, isSecondaryAlignment = NA, isDuplicate = FALSE,
#    isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE
#  )
#  scan_param <- ScanBamParam(
#    flag = scanflag,
#    what = c("qname", "flag", "cigar"),
#    tag = "NM"
#  )
#  aln <- scanBam(bamf, param = scan_param)
#  aln_dt <- data.table(
#    qname = aln[[1]]$qname,
#    cigar = aln[[1]]$cigar,
#    flag = aln[[1]]$flag,
#    nm = unlist(aln[[1]]$tag)
#  )
#  aln_dt[, read_idx := ifelse(
#    bamFlagAsBitMatrix(as.integer(flag))[7] == 1, 1, 2
#  ),
#  by = seq_len(nrow(aln_dt))
#  ]
#  if (nrow(aln_dt) == 0) {
#    print(paste(
#      "[ERROR] Found no alignments in the give BAM: ", bam,
#      sep = ""
#    ))
#    quit(status = 1)
#  }
#  cigar <- NULL # this is to avoid "no visible binding for cigar"
#  n_ins <- n_del <- n_mm <- n_ecnt <- NULL
#  nread_per_frag <- NULL
#  aln_dt[, ":="(
#    n_ins = length(
#      grep(pattern = "I", unlist(strsplit(cigar, "")))
#    ),
#    n_del = length(
#      grep(pattern = "D", unlist(strsplit(cigar, "")))
#    )
#  ), by = seq_len(nrow(aln_dt))]
#  aln_dt[, "n_mm" := nm - apply(.SD, 1, get_indel_length), .SDcols = "cigar"]
#  aln_dt[, "n_ecnt" := n_mm + n_ins + n_del]
#  aln_dt <- aln_dt[n_ecnt <= min_necnt]
#  aln_dt[, "nread_per_frag" := .N, by = "qname"]
#  # because we are including secondary alignments, here using mod operator
#  aln_dt <- aln_dt[nread_per_frag %% 2 == 0]
#  filter <- S4Vectors::FilterRules(
#    list(function(x) x$qname %in% aln_dt$qname)
#  )
#  obam <- filterBam(bamf, obam, filter = filter, param = scan_param)
# }

# init_loh_report <- function(a1, a2) {
#  list(
#    "HLA_A1" = a1,
#    "HLA_A2" = a2,
#    "HLA_A1_CN" = NaN,
#    "HLA_A1_CN_Lower" = NaN,
#    "HLA_A1_CN_Upper" = NaN,
#    "HLA_A2_CN" = NaN,
#    "HLA_A2_CN_Lower" = NaN,
#    "HLA_A2_CN_Upper" = NaN,
#    "HLA_A1_Median_LogR" = NaN,
#    "HLA_A2_Median_LogR" = NaN,
#    "HLA_A1_MM_Median_LogR" = NaN,
#    "HLA_A2_MM_Median_LogR" = NaN,
#    "MM_LogR_Paired_Pvalue" = NaN,
#    "Median_BAF" = NaN,
#    "Num_MM" = as.integer(0),
#    "Num_Bins" = as.integer(0),
#    "Num_Bins_Used" = as.integer(0),
#    "Num_CN_Loss_Supporting_Bins" = as.integer(0)
#  )
# }

# extract_allele_coverage <- function(allele, bam, hlaref, min_dp = 0) {
#  dt <- fread(text = system2(
#    command = "samtools",
#    args = c(
#      "mpileup", "-f", hlaref,
#      "--rf", 2, "--ff", 3584, "-Q", 20, "-r", allele, bam
#    ),
#    stdout = TRUE,
#    stderr = FALSE,
#    wait = TRUE
#  ), select = c(1, 2, 3, 4))
#  names(dt) <- c("seqnames", "pos", "nucleotide", "count")
#  dt <- dt[count > min_dp]
#  dt
# }

# get_allele_coverage <- function(allele, bam, min_dp = 0) {
#  print(paste(
#    "[INFO] Get coverage for ", allele, " from ", bam,
#    sep = ""
#  ))
#  seqinfo <- extract_seqinfo_from_bam(bam = bam)
#  allele_seq_ln <- seqinfo[which(names(seqinfo) == allele)]
#  allele_to_scan <- GenomicRanges::GRanges(
#    seqnames = allele,
#    ranges = IRanges::IRanges(start = 1, end = allele_seq_ln)
#  )
#  scanflag <- scanBamFlag(
#    isProperPair = TRUE, isSecondaryAlignment = FALSE, isDuplicate = FALSE,
#    isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE
#  )
#  scan_param <- ScanBamParam(
#    flag = scanflag,
#    mapqFilter = 20,
#    which = allele_to_scan,
#  )
#  pileup_param <- PileupParam(
#    min_mapq = 20, distinguish_strands = FALSE, max_depth = 9999
#  )
#  p_dt <- setDT(
#    pileup(
#      file = bam, scanBamParam = scan_param, pileupParam = pileup_param
#    )
#  )
#  # pileup return seqnames as factor
#  p_dt <- p_dt[count > min_dp]
#  p_dt[, seqnames := as.character(seqnames)]
#  p_dt[, which_label := NULL]
#  p_dt
# }

# prep_allelic_cov <- function(t_dt, n_dt, bin_dt, multfactor) {
#  cov_dt <- combine_tn_cov(t_dt = t_dt, n_dt = n_dt)
#  cov_dt <- bin_allele_cov(cov_dt = cov_dt, bin_dt = bin_dt)
#  cov_dt <- estimate_logr(cov_dt = cov_dt, multfactor = multfactor)
#  cov_dt
# }

# combine_tn_cov <- function(t_dt, n_dt) {
#  a_t <- unique(t_dt$seqnames)
#  a_n <- unique(n_dt$seqnames)
#  if (length(a_t) != 1 || length(a_n) != 1) {
#    print("[ERROR] Only one seqname is expected in tumor and normal tables")
#    quit(status = 1)
#  }
#  if (a_t != a_n) {
#    print("[ERROR] To be combined tumor and normal tables have diff seqnames")
#    quit(status = 1)
#  }
#  setkey(t_dt, seqnames, pos, nucleotide)
#  setkey(n_dt, seqnames, pos, nucleotide)
#  allele_cov_dt <- t_dt[n_dt]
#  allele_cov_dt[, ":="(t_dp = count, n_dp = i.count)]
#  allele_cov_dt[, ":="(count = NULL, i.count = NULL)]
#
#  allele_cov_dt
# }

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

# make_bins <- function(aln, bin_size = 150) {
#  aln$aln_dt[, bin := (pa_pos - 1) %/% bin_size]
#  aln$aln_dt[, bin := bin + 1, by = bin]
#  if (nrow(aln$aln_dt[bin == max(bin)]) < bin_size) {
#    aln$aln_dt[bin == max(bin), bin := max(bin) - 1]
#  }
#
#  aln$aln_dt[, ":="(
#    a1_start = min(a1_pos),
#    a1_end = max(a1_pos),
#    a2_start = min(a2_pos),
#    a2_end = max(a2_pos)
#  ), by = bin]
#  bin_dt <- unique(aln$aln_dt, by = "bin")
#
#  list(
#    a1_bin = bin_dt[, c("a1_seqnames", "a1_start", "a1_end", "bin")],
#    a2_bin = bin_dt[, c("a2_seqnames", "a2_start", "a2_end", "bin")]
#  )
# }

# bin_allele_cov <- function(cov_dt, bin_dt) {
#  names(bin_dt) <- gsub("^.*_", "", names(bin_dt))
#  # make granges for allele coverage dt and bins
#  tmp_gr <- GenomicRanges::makeGRangesFromDataFrame(
#    cov_dt,
#    start.field = "pos", end.field = "pos"
#  )
#  bin_gr <- GenomicRanges::makeGRangesFromDataFrame(
#    bin_dt,
#    keep.extra.columns = TRUE
#  )
#  ovl <- GenomicRanges::findOverlaps(tmp_gr, bin_gr)
#  cov_dt[S4Vectors::queryHits(ovl), "bin"] <- bin_dt[
#    S4Vectors::subjectHits(ovl), "bin"
#  ]
#  cov_dt[, ":="(
#    bin_t_dp = as.numeric(median(t_dp, na.rm = TRUE)),
#    bin_n_dp = as.numeric(median(n_dp, na.rm = TRUE))
#  ), by = "bin"]
#  cov_dt
# }

# estimate_logr <- function(cov_dt, multfactor) {
#  cov_dt[, logR := log2(t_dp / n_dp * multfactor)]
#  # FIXME: i can also use bin_t_dp and bin_n_dp calculate
#  # two should be similar
#  cov_dt[, bin_logR := median(logR, na.rm = TRUE), by = "bin"]
#  cov_dt
# }

# estimate_binned_logr <- function(a1_dt, a2_dt, multfactor) {
#  bin_dt <- merge(
#    unique(a1_dt, by = "a1_bin"),
#    unique(a2_dt, by = "a2_bin"),
#    by.x = c("a1_bin"), by.y = c("a2_bin"),
#    all.x = TRUE, all.y = TRUE
#  )
#  binned_names <- names(bin_dt)[which(grepl("bin", names(bin_dt)))]
#  bin_dt <- bin_dt[, ..binned_names]
#  bin_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
#  bin_dt[, logR_combined_bin := log2(
#    (a1_bin_t_dp + a2_bin_t_dp) /
#      (a1_bin_n_dp + a2_bin_n_dp) * multfactor
#  )]
#  bin_dt
# }

# prep_mm_cov <- function(mm, a1_dt, a2_dt) {
#  mm_est_dt <- data.table(a1_pos = mm$diffSeq1, a2_pos = mm$diffSeq2)
#  mm_est_dt <- mm_est_dt[
#    a1_pos %in% a1_dt$a1_pos & a2_pos %in% a2_dt$a2_pos
#  ]
#  if (nrow(mm_est_dt) == 0) {
#    return(NULL)
#  }
#  simply_cols <- names(a1_dt)[which(!grepl("bin", names(a1_dt)))]
#  simply_cols <- c(simply_cols, "a1_bin")
#  mm_est_dt <- merge(
#    mm_est_dt, a1_dt[, ..simply_cols],
#    by = "a1_pos", all.x = TRUE
#  )
#  mm_est_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
#  if (nrow(mm_est_dt) == 0) {
#    return(NULL)
#  }
#  simply_cols <- names(a2_dt)[which(!grepl("bin", names(a2_dt)))]
#  mm_est_dt <- merge(
#    mm_est_dt, a2_dt[, ..simply_cols],
#    by = "a2_pos", all.x = TRUE
#  )
#  if (nrow(mm_est_dt) == 0) {
#    return(NULL)
#  }
#
#  mm_est_dt
# }

# estimate_baf <- function(mm_dt, capture_bias) {
#  mm_dt[, baf := a1_t_dp / (a1_t_dp + a2_t_dp)] # nolint
#  mm_dt[, baf_correct := baf / capture_bias]
#  n_odd_baf <- nrow(mm_dt[baf_correct > 1])
#  if (n_odd_baf > 0) {
#    print("[WARN] Found mismatch sites with corrected BAF larger than 1")
#    print(mm_dt[baf_correct > 1])
#    print("[WARN] Please check if two alleles have consistent coverage")
#    print("[WARM] Force setting BAF to 1")
#  }
#  mm_dt[, baf_correct := ifelse(baf_correct > 1, 1, baf_correct)]
#  mm_dt
# }

# estimate_cn <- function(mm_dt, bin_dt) {
#  req_cols <- c("bin", "logR_combined_bin")
#  miss_cols <- req_cols[which(!req_cols %in% names(bin_dt))]
#  if (length(miss_cols) > 0) {
#    print(paste(
#      "[ERROR] Miss ",
#      paste(miss_cols, collapse = ","),
#      " columns in bin_dt",
#      sep = ""
#    ))
#    print("[ERROR] Cannot continue estimate BAF")
#    quit(status = 1)
#  }
#  mm_dt <- merge(
#    mm_dt,
#    bin_dt[, c("bin", "logR_combined_bin")],
#    by = "bin",
#    all.x = TRUE
#  )
#  mm_dt[, "a1_cn" :=
#    (purity - 1 + baf_correct * 2^(logR_combined_bin / gamma) * # nolint
#      ((1 - purity) * 2 + purity * ploidy)) / purity, ]
#  mm_dt[, "a2_cn" :=
#    (purity - 1 - (baf_correct - 1) * 2^(logR_combined_bin / gamma) * # nolint
#      ((1 - purity) * 2 + purity * ploidy)) / purity, ]
#  mm_dt[baf_correct > 1, ":="(a1_cn = NaN, a2_cn = NaN)]
#  mm_dt
# }

# estimate_cn_conf <- function(cn_dt, which) {
#  col <- NULL
#  pattern <- paste(which, "bin_cn", sep = "_")
#  col <- names(cn_dt)[which(grepl(pattern, names(cn_dt)))]
#  if (is.null(col) || length(col) == 0) {
#    print(paste(
#      "[ERROR] Failed to find cn column in cn_dt for ",
#      which, " allele",
#      sep = ""
#    ))
#    quit(status = 1)
#  }
#  cn_test <- t_test_with_na(cn_dt[[col]])
#  cn_est_conf <- cn_test$conf.int
#  cn_est_lower <- cn_est_conf[1]
#  cn_est_upper <- cn_est_conf[2]
#  list(est_lower = cn_est_lower, est_upper = cn_est_upper)
# }

# dump_intermedia_tables <- function(
#    out, dp_info, mm,
#    cov_dt, bin_dt = NULL, mm_dt = NULL) {
#  to_save <- list(
#    mm = mm,
#    cov_dt = cov_dt,
#    bin_dt = bin_dt,
#    mm_dt = mm_dt,
#    dp_info = dp_info
#  )
#  saveRDS(to_save, out)
# }

# estimate_dp <- function(bam, alleles) {
#  count_n_reads_from_bam(
#    bam = bam, which = alleles, tagfilter = list(NM = c(0, 1))
#  )
# }

# estimate_capture_bias <- function(n_a1_cov, n_a2_cov, aln) {
#  n_cov_dt <- merge(n_a1_cov, aln, by.x = "pos", by.y = "a1_pos")
#  n_cov_dt <- merge(n_a2_cov, n_cov_dt, by.x = "pos", by.y = "a2_pos")
#  n_cov_dt[, cov_ratio := count.x / count.y]
#  median(n_cov_dt$cov_ratio, na.rm = TRUE)
# }

# call_hla_loh <- function(
#    dt, tbam, nbam, hlaref, outdir,
#    purity, ploidy, min_dp, min_necnt,
#    tid = "example_tumor", nid = "example_normal", gamma = 1) {
#  a1 <- dt$A1
#  a2 <- dt$A2
#  alleles_str <- paste(c(a1, a2), collapse = " and ")
#  print(paste("[INFO] Analyze LOH for ", alleles_str, sep = ""))
#
#  report <- init_loh_report(a1, a2)
#
#  n_seq_depth <- estimate_dp(bam = nbam, alleles = c(a1, a2))
#  t_seq_depth <- estimate_dp(bam = tbam, alleles = c(a1, a2))
#  dp_info <- list(n_seq_depth = n_seq_depth, t_seq_depth = t_seq_depth)
#  multfactor <- n_seq_depth / t_seq_depth
#  print(paste(
#    "[INFO] Align sequences between ", alleles_str,
#    sep = ""
#  ))
#  mm <- get_mismatches_bw_alleles(a1, a2, hlaref)
#  if (is.null(mm)) {
#    print("[INFO] no call will be made. Move to the next HLA gene")
#    return(report)
#  }
#  report$Num_MM <- length(mm$diffSeq1)
#
#  print("[INFO] Make bins accounting for alignment start position")
#  bins <- make_bins(aln = mm)
#  report$Num_Bins <- nrow(bins$a1_bin)
#
#  t_a1_cov <- extract_allele_coverage(allele = a1, bam = tbam, hlaref = hlaref)
#  t_a2_cov <- extract_allele_coverage(allele = a2, bam = tbam, hlaref = hlaref)
#  n_a1_cov <- extract_allele_coverage(
#    allele = a1, bam = nbam, hlaref = hlaref, min_dp = min_dp
#  )
#  n_a2_cov <- extract_allele_coverage(
#    allele = a2, bam = nbam, hlaref = hlaref, min_dp = min_dp
#  )
#  t_a1_cov <- t_a1_cov[pos %in% n_a1_cov$pos]
#  t_a2_cov <- t_a2_cov[pos %in% n_a2_cov$pos]
#  capture_bias <- estimate_capture_bias(
#    n_a1_cov = n_a1_cov, n_a2_cov = n_a2_cov, aln = mm$aln_dt
#  )
#  if (nrow(n_a1_cov) == 0 || nrow(n_a2_cov) == 0) {
#    print("[INFO] Found no coverage for either allele in normal")
#    print("[INFO] Move to next HLA gene")
#    return(report)
#  }
#  print(paste(
#    "[INFO] Prepare coverage table for allele ", a1,
#    sep = ""
#  ))
#  a1_cov_dt <- prep_allelic_cov(
#    t_dt = t_a1_cov, n_dt = n_a1_cov,
#    bin_dt = bins$a1_bin, multfactor = multfactor
#  )
#  names(a1_cov_dt) <- paste("a1_", names(a1_cov_dt), sep = "")
#  print(paste(
#    "[INFO] Prepare coverage table for allele ", a2,
#    sep = ""
#  ))
#  a2_cov_dt <- prep_allelic_cov(
#    t_dt = t_a2_cov, n_dt = n_a2_cov,
#    bin_dt = bins$a2_bin, multfactor = multfactor
#  )
#  names(a2_cov_dt) <- paste("a2_", names(a2_cov_dt), sep = "")
#
#  bin_dt <- estimate_binned_logr(
#    a1_dt = a1_cov_dt, a2_dt = a2_cov_dt, multfactor = multfactor
#  )
#  bin_dt[, cn_loss_test_bin := apply(
#    .SD, 1, test_cn_loss
#  ),
#  .SDcols = c("a1_bin_t_dp", "a1_bin_n_dp", "a2_bin_t_dp", "a2_bin_n_dp")
#  ]
#  setkey(bin_dt, bin)
#
#  mm_dt <- prep_mm_cov(mm = mm, a1_dt = a1_cov_dt, a2_dt = a2_cov_dt)
#  if (is.null(mm_dt)) {
#    print("[INFO] No coverage info found for both alleles at mm sites")
#    print("[INFO] Cannot estimate copy numbers. Moving to next HLA gene")
#    return(report)
#  }
#  mm_dt <- estimate_baf(mm_dt = mm_dt, capture_bias = capture_bias)
#  paired_t_test <- t.test(mm_dt$a1_logR, mm_dt$a2_logR, paired = TRUE)
#  report$MM_LogR_Paired_Pvalue <- paired_t_test$p.value
#  report$Median_BAF <- median(mm_dt$baf_correct, na.rm = TRUE)
#
#  print("[INFO] Estimate copy number at mismatch sites")
#  mm_dt <- estimate_cn(mm_dt = mm_dt, bin_dt = bin_dt)
#  mm_dt[, ":="(
#    a1_bin_cn = median(a1_cn, na.rm = TRUE),
#    a2_bin_cn = median(a2_cn, na.rm = TRUE)
#  ),
#  by = "bin"
#  ]
#  # FIXME: I should get this right after finished making bins
#  bins_no_mm <- bin_dt$bin[which(!bin_dt$bin %in% unique(mm_dt$bin))]
#  report$Num_Bins_Used <- report$Num_Bins - length(bins_no_mm)
#  report$HLA_A1_Median_LogR <- median(
#    a1_cov_dt[!a1_bin %in% bins_no_mm]$a1_logR,
#    na.rm = TRUE
#  )
#  report$HLA_A2_Median_LogR <- median(
#    a2_cov_dt[!a2_bin %in% bins_no_mm]$a2_logR,
#    na.rm = TRUE
#  )
#  report$HLA_A1_MM_Median_LogR <- median(mm_dt$a1_logR, na.rm = TRUE)
#  report$HLA_A2_MM_Median_LogR <- median(mm_dt$a2_logR, na.rm = TRUE)
#
#  cn_dt <- unique(mm_dt, by = "bin")
#  setkey(cn_dt, bin)
#  report$HLA_A1_CN <- round(median(cn_dt$a1_bin_cn, na.rm = TRUE), 4)
#  report$HLA_A2_CN <- round(median(cn_dt$a2_bin_cn, na.rm = TRUE), 4)
#
#  cn_est_conf <- estimate_cn_conf(cn_dt = cn_dt, which = "a1")
#  report$HLA_A1_CN_Lower <- round(cn_est_conf$est_lower, 4)
#  report$HLA_A1_CN_Upper <- round(cn_est_conf$est_upper, 4)
#  cn_est_conf <- estimate_cn_conf(cn_dt = cn_dt, which = "a2")
#  report$HLA_A2_CN_Lower <- round(cn_est_conf$est_lower, 4)
#  report$HLA_A2_CN_Upper <- round(cn_est_conf$est_upper, 4)
#  bin_dt <- cn_dt[, c("bin", "a1_bin_cn", "a2_bin_cn")][bin_dt]
#  report$Num_CN_Loss_Supporting_Bins <- nrow(
#    bin_dt[cn_loss_test_bin <= 0.01 & !bin %in% bins_no_mm]
#  )
#  rm(cn_dt)
#
#  hla_gene <- gsub("(_|\\*)*(([0-9])+(_|:)*)+", "", a1)
#  hla_gene <- tolower(hla_gene)
#  out_rds <- file.path(outdir, paste(hla_gene, ".rds", sep = ""))
#  names(a1_cov_dt) <- gsub("a1_", "", names(a1_cov_dt))
#  names(a2_cov_dt) <- gsub("a2_", "", names(a2_cov_dt))
#  dump_intermedia_tables(
#    out = out_rds,
#    dp_info = dp_info,
#    mm = mm,
#    cov_dt = rbind(a1_cov_dt, a2_cov_dt),
#    bin_dt = bin_dt,
#    mm_dt = mm_dt
#  )
#
#  report
# }

initialize_libs <- function() {
  pkg_name <- "hlalohReforged"
  libpaths <- file.path(.libPaths(), pkg_name)
  lib_pattern <- "(bamer|cli|loh|pairwise_aln|pathio).R"
  rscripts <- list.files(libpaths, lib_pattern, full.names = TRUE, recursive = TRUE)
  if (length(rscripts) == 0) {
    print("[ERROR] Found no related R script libraries to run hlalohReforged")
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
# gamma <- 1
#
# args <- parse_cmd()
#
# parse_file_path(file = args$tbam)
# parse_file_path(file = args$nbam)
# parse_dir_path(dir = args$outdir, create = TRUE)
#
# tid <- extract_rgsm_from_bam_header(bam = args$tbam)
# nid <- extract_rgsm_from_bam_header(bam = args$nbam)
#
# purity <- 0.5
# ploidy <- 2
# if (!is.null(args$tstates)) {
#  tstates <- extract_tstates(tstate_est_file = args$tstates)
#  if (!is.na(tstates$ploidy) && !is.na(tstates$purity)) {
#    ploidy <- tstates$ploidy
#    purity <- tstates$purity
#  }
# }
# print(paste("[INFO] Purity = ", purity, " Ploidy = ", ploidy, sep = ""))
#
# alleles_n <- extract_seqinfo_from_bam(bam = args$nbam)
# alleles_t <- extract_seqinfo_from_bam(bam = args$tbam)
## extract function returns named integer vector
## here we need names of that vector to get allele name
# alleles_n <- sort(names(alleles_n))
# alleles_t <- sort(names(alleles_t))
# if (!all.equal(alleles_n, alleles_t)) {
#  print("[ERROR] Different set of alleles detected in normal and tumor BAMs")
#  print("[ERROR] Alleles in normal BAM: ", paste(alleles_n, collapse = "|"))
#  print("[ERROR] Alleles in tumor BAM: ", paste(alleles_t, collapse = "|"))
#  quit(status = 1)
# }
#
# alleles_dt <- data.table(Alleles = alleles_n)
# alleles_dt[, HLAGene := tstrsplit(Alleles, "_", keep = 2)]
# alleles_dt[, HLAGene := paste("hla_", HLAGene, sep = "")]
# alleles_dt[, Alleles := paste(Alleles, collapse = ","), by = HLAGene]
# alleles_dt[, c("A1", "A2") := tstrsplit(Alleles, ","), by = HLAGene]
# alleles_dt[, Alleles := NULL]
# alleles_dt <- unique(alleles_dt, by = "HLAGene")
#
## filter input bams by ecnt
# filt_nbam <- file.path(args$outdir, paste(nid, ".filt.bam", sep = ""))
# filt_tbam <- file.path(args$outdir, paste(tid, ".filt.bam", sep = ""))
# if (!file.exists(filt_nbam) || !file.exists(filt_tbam)) {
#  filter_bam_by_ecnt(
#    bam = args$nbam, obam = filt_nbam, min_necnt = args$min_necnt
#  )
#  filter_bam_by_ecnt(
#    bam = args$tbam, obam = filt_tbam, min_necnt = args$min_necnt
#  )
# }
#
# loh_res_dt <- alleles_dt[, call_hla_loh(
#  .SD,
#  tbam = filt_tbam, nbam = filt_nbam, hlaref = args$hlaref,
#  outdir = args$outdir, purity = purity, ploidy = ploidy,
#  min_dp = args$min_cov, min_necnt = args$min_nm,
#  tid = tid, nid = nid
# ), by = "HLAGene"]
# out_res <- file.path(args$outdir, paste(args$subject, ".loh.res.tsv", sep = ""))
# fwrite(loh_res_dt, out_res, sep = "\t", row.names = FALSE, quote = FALSE)
