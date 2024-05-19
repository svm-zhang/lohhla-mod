require(data.table)

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

bin_allele_cov <- function(cov_dt, bin_dt) {
  names(bin_dt) <- gsub("^.*_", "", names(bin_dt))
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

make_bins <- function(aln, bin_size = 150) {
  aln$aln_dt[, bin := (pa_pos - 1) %/% bin_size]
  aln$aln_dt[, bin := bin + 1, by = bin]
  if (nrow(aln$aln_dt[bin == max(bin)]) < bin_size) {
    aln$aln_dt[bin == max(bin), bin := max(bin) - 1]
  }

  aln$aln_dt[, ":="(
    a1_start = min(a1_pos),
    a1_end = max(a1_pos),
    a2_start = min(a2_pos),
    a2_end = max(a2_pos)
  ), by = bin]
  bin_dt <- unique(aln$aln_dt, by = "bin")

  list(
    a1_bin = bin_dt[, c("a1_seqnames", "a1_start", "a1_end", "bin")],
    a2_bin = bin_dt[, c("a2_seqnames", "a2_start", "a2_end", "bin")]
  )
}

estimate_capture_bias <- function(n_a1_cov, n_a2_cov, aln) {
  n_cov_dt <- merge(n_a1_cov, aln, by.x = "pos", by.y = "a1_pos")
  n_cov_dt <- merge(n_a2_cov, n_cov_dt, by.x = "pos", by.y = "a2_pos")
  n_cov_dt[, cov_ratio := count.x / count.y]
  median(n_cov_dt$cov_ratio, na.rm = TRUE)
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

estimate_baf <- function(mm_dt, capture_bias) {
  mm_dt[, baf := a1_t_dp / (a1_t_dp + a2_t_dp)] # nolint
  mm_dt[, baf_correct := baf / capture_bias]
  n_odd_baf <- nrow(mm_dt[baf_correct > 1])
  if (n_odd_baf > 0) {
    print("[WARN] Found mismatch sites with corrected BAF larger than 1")
    print(mm_dt[baf_correct > 1])
    print("[WARN] Please check if two alleles have consistent coverage")
    print("[WARM] Force setting BAF to 1")
  }
  mm_dt[, baf_correct := ifelse(baf_correct > 1, 1, baf_correct)]
  mm_dt
}

estimate_cn <- function(mm_dt, bin_dt, purity, ploidy, gamma = 1) {
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
    bin_dt[, c("bin", "logR_combined_bin")],
    by = "bin",
    all.x = TRUE
  )
  mm_dt[, "a1_cn" :=
    (purity - 1 + baf_correct * 2^(logR_combined_bin / gamma) * # nolint
      ((1 - purity) * 2 + purity * ploidy)) / purity, ]
  mm_dt[, "a2_cn" :=
    (purity - 1 - (baf_correct - 1) * 2^(logR_combined_bin / gamma) * # nolint
      ((1 - purity) * 2 + purity * ploidy)) / purity, ]
  mm_dt[baf_correct > 1, ":="(a1_cn = NaN, a2_cn = NaN)]
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

estimate_dp <- function(bam, alleles) {
  count_n_reads_from_bam(
    bam = bam, which = alleles, tagfilter = list(NM = c(0, 1))
  )
}

dump_intermedia_tables <- function(
    out, dp_info, mm,
    cov_dt, bin_dt = NULL, mm_dt = NULL) {
  to_save <- list(
    mm = mm,
    cov_dt = cov_dt,
    bin_dt = bin_dt,
    mm_dt = mm_dt,
    dp_info = dp_info
  )
  saveRDS(to_save, out)
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
    "Num_Bins_Used" = as.integer(0),
    "Num_CN_Loss_Supporting_Bins" = as.integer(0)
  )
}

call_hla_loh <- function(
    dt, tbam, nbam, hlaref, outdir,
    purity, ploidy, min_dp, min_necnt, gamma = 1) {
  a1 <- dt$A1
  a2 <- dt$A2
  alleles_str <- paste(c(a1, a2), collapse = " and ")
  print(paste("[INFO] Analyze LOH for ", alleles_str, sep = ""))

  report <- init_loh_report(a1, a2)

  n_seq_depth <- estimate_dp(bam = nbam, alleles = c(a1, a2))
  t_seq_depth <- estimate_dp(bam = tbam, alleles = c(a1, a2))
  dp_info <- list(n_seq_depth = n_seq_depth, t_seq_depth = t_seq_depth)
  multfactor <- n_seq_depth / t_seq_depth
  print(paste(
    "[INFO] Align sequences between ", alleles_str,
    sep = ""
  ))
  mm <- get_mismatches_bw_alleles(a1, a2, hlaref)
  if (is.null(mm)) {
    print("[INFO] no call will be made. Move to the next HLA gene")
    return(report)
  }
  report$Num_MM <- length(mm$diffSeq1)

  print("[INFO] Make bins accounting for alignment start position")
  bins <- make_bins(aln = mm)
  report$Num_Bins <- nrow(bins$a1_bin)

  t_a1_cov <- extract_allele_coverage(allele = a1, bam = tbam, hlaref = hlaref)
  t_a2_cov <- extract_allele_coverage(allele = a2, bam = tbam, hlaref = hlaref)
  n_a1_cov <- extract_allele_coverage(
    allele = a1, bam = nbam, hlaref = hlaref, min_dp = min_dp
  )
  n_a2_cov <- extract_allele_coverage(
    allele = a2, bam = nbam, hlaref = hlaref, min_dp = min_dp
  )
  t_a1_cov <- t_a1_cov[pos %in% n_a1_cov$pos]
  t_a2_cov <- t_a2_cov[pos %in% n_a2_cov$pos]
  capture_bias <- estimate_capture_bias(
    n_a1_cov = n_a1_cov, n_a2_cov = n_a2_cov, aln = mm$aln_dt
  )
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
    bin_dt = bins$a1_bin, multfactor = multfactor
  )
  names(a1_cov_dt) <- paste("a1_", names(a1_cov_dt), sep = "")
  print(paste(
    "[INFO] Prepare coverage table for allele ", a2,
    sep = ""
  ))
  a2_cov_dt <- prep_allelic_cov(
    t_dt = t_a2_cov, n_dt = n_a2_cov,
    bin_dt = bins$a2_bin, multfactor = multfactor
  )
  names(a2_cov_dt) <- paste("a2_", names(a2_cov_dt), sep = "")

  bin_dt <- estimate_binned_logr(
    a1_dt = a1_cov_dt, a2_dt = a2_cov_dt, multfactor = multfactor
  )
  bin_dt[, cn_loss_test_bin := apply(
    .SD, 1, test_cn_loss
  ),
  .SDcols = c("a1_bin_t_dp", "a1_bin_n_dp", "a2_bin_t_dp", "a2_bin_n_dp")
  ]
  setkey(bin_dt, bin)

  mm_dt <- prep_mm_cov(mm = mm, a1_dt = a1_cov_dt, a2_dt = a2_cov_dt)
  if (is.null(mm_dt)) {
    print("[INFO] No coverage info found for both alleles at mm sites")
    print("[INFO] Cannot estimate copy numbers. Moving to next HLA gene")
    return(report)
  }
  mm_dt <- estimate_baf(mm_dt = mm_dt, capture_bias = capture_bias)
  paired_t_test <- t.test(mm_dt$a1_logR, mm_dt$a2_logR, paired = TRUE)
  report$MM_LogR_Paired_Pvalue <- paired_t_test$p.value
  report$Median_BAF <- median(mm_dt$baf_correct, na.rm = TRUE)

  print("[INFO] Estimate copy number at mismatch sites")
  mm_dt <- estimate_cn(
    mm_dt = mm_dt, bin_dt = bin_dt,
    purity = purity, ploidy = ploidy, gamma = gamma
  )
  mm_dt[, ":="(
    a1_bin_cn = median(a1_cn, na.rm = TRUE),
    a2_bin_cn = median(a2_cn, na.rm = TRUE)
  ),
  by = "bin"
  ]
  # FIXME: I should get this right after finished making bins
  bins_no_mm <- bin_dt$bin[which(!bin_dt$bin %in% unique(mm_dt$bin))]
  report$Num_Bins_Used <- report$Num_Bins - length(bins_no_mm)
  report$HLA_A1_Median_LogR <- median(
    a1_cov_dt[!a1_bin %in% bins_no_mm]$a1_logR,
    na.rm = TRUE
  )
  report$HLA_A2_Median_LogR <- median(
    a2_cov_dt[!a2_bin %in% bins_no_mm]$a2_logR,
    na.rm = TRUE
  )
  report$HLA_A1_MM_Median_LogR <- median(mm_dt$a1_logR, na.rm = TRUE)
  report$HLA_A2_MM_Median_LogR <- median(mm_dt$a2_logR, na.rm = TRUE)

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
  report$Num_CN_Loss_Supporting_Bins <- nrow(
    bin_dt[cn_loss_test_bin <= 0.01 & !bin %in% bins_no_mm]
  )
  rm(cn_dt)

  hla_gene <- gsub("(_|\\*)*(([0-9])+(_|:)*)+", "", a1)
  hla_gene <- tolower(hla_gene)
  out_rds <- file.path(outdir, paste(hla_gene, ".rds", sep = ""))
  names(a1_cov_dt) <- gsub("a1_", "", names(a1_cov_dt))
  names(a2_cov_dt) <- gsub("a2_", "", names(a2_cov_dt))
  dump_intermedia_tables(
    out = out_rds,
    dp_info = dp_info,
    mm = mm,
    cov_dt = rbind(a1_cov_dt, a2_cov_dt),
    bin_dt = bin_dt,
    mm_dt = mm_dt
  )

  report
}
