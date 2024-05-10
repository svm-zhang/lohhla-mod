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
    metavar = "FILE", type = "character", required = TRUE,
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

test_cn_loss_using_log_odds_ratio <- function(x) {
  if(anyNA(x)) {
    return(NA)
  }
  x <- matrix(round(x), nrow=2)
  test <- fisher.test(x)

  test$p.value
}

normalize_hla_seqname <- function(dt) {

  if (! "seqnames" %in% names(dt)) {
    print("[ERROR] Cannot find seqnames column in the given data table")
    quit(status = 1)
  }
  dt[, seqnames := gsub("hla_", "", seqnames)]
  dt[, hla_gene := unlist(strsplit(seqnames, "_"))[1], by=seq_len(nrow(dt))]
  dt[, seqnames := gsub(
    paste(hla_gene, "_", sep=""),
    paste(toupper(hla_gene), "*", sep=""),
    seqnames
  ),
    by = seq_len(nrow(dt))
  ]
  dt[, seqnames := gsub("_", ":", seqnames)]
  dt[, hla_gene := NULL]

  dt
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

get_read_group_from_bam <- function(bam) {
  bf <- BamFile(bam)

  rg <- scanBamHeader(bf)$text$`@RG`
  if (is.null(rg) || is.na(rg)) {
    print(paste("[ERROR] Cannot find the read group (@RG) in the BAM file: ",
      bam,
      sep = ""
    ))
    quit(status = 1)
  }
  rg
}

get_sm_from_rg <- function(rg) {
  sm <- rg[which(grepl("^SM", rg))]
  if (length(sm) == 0) {
    print(paste("[ERROR] Cannot find the SM in the @RG: ",
      paste(rg, collapse = "\t"),
      sep = ""
    ))
    quit(status = 1)
  }
  sm <- gsub("^SM:", "", sm)
  sm
}

get_hla_alleles <- function(hla_res, hla_ref_fasta) {
  print("[INFO] Loading HLA typing result")
  expected_cols_in_hla_res <- c(
    "SampleID", "HLA-A_1", "HLA-A_2",
    "HLA-B_1", "HLA-B_2", "HLA-C_1", "HLA-C_2"
  )
  hla_res_dt <- fread(hla_res)
  print("[INFO] Loading HLA typing result [DONE]")

  # to make sure any missing data to be replaced by dot
  # you can trust anything, sorry
  for (i in names(hla_res_dt)) {
    hla_res_dt[is.na(get(i)), (i) := "."]
  }

  # make sure we have the right columns to work with
  if (!all.equal(names(hla_res_dt), expected_cols_in_hla_res)) {
    print("[ERROR] Missing expected columns in the HLA typing result table")
    print(paste(
      "[ERROR] Expected columns are: ", expected_cols_in_hla_res,
      sep = ""
    ))
    quit(status = 1)
  }

  hla_res_dt <- melt(
    data = hla_res_dt,
    id.vars = "SampleID",
    variable.name = "HLA_Alleles", value.name = "Allele_Type"
  )
  hla_res_dt[, "HLA_Alleles" := ifelse(grepl("^HLA-A", HLA_Alleles), "HLA-A", # nolint
    ifelse(grepl("^HLA-B", HLA_Alleles), "HLA-B",
      ifelse(grepl("^HLA-C", HLA_Alleles), "HLA-C", "Other")
    )
  )]
  hla_res_dt <- hla_res_dt[,
    .(Allele_Type = .(Allele_Type)), # nolint
    by = .(HLA_Alleles) # nolint
  ]
  hla_alleles <- lapply(
    split(hla_res_dt, as.factor(hla_res_dt$HLA_Alleles)),
    function(x) unlist(x$Allele_Type)
  )

  for (hla_gene in names(hla_alleles)) {
    hla_alleles[[hla_gene]] <- hla_alleles[[hla_gene]][
      !is.na(hla_alleles[[hla_gene]])
    ]
    alleles <- unique(hla_alleles[[hla_gene]])
    if (length(alleles) == 1) {
      if (alleles == ".") {
        print(paste(
          "[INFO] ", hla_gene, " allele missing from the HLA typing result",
          sep = ""
        ))
      } else {
        print(paste(
          "[INFO] ", hla_gene, " is homozygous for ", alleles,
          sep = ""
        ))
      }
      hla_alleles[[hla_gene]] <- NULL
      print(paste("[INFO] ", hla_gene, " wont be analyzed for LOH", sep = ""))
    } else if (length(alleles) > 2) {
      print(paste(
        "[ERROR] Only expect diploid ", hla_gene,
        " gene. But see ", alleles,
        sep = ""
      ))
      quit(status = 1)
    } else {
      hla_alleles[[hla_gene]] <- check_if_alleles_in_ref(
        alleles = alleles, hla_ref_fasta = hla_ref_fasta
      )
    }
  }

  if (length(hla_alleles) == 0) {
    print("[INFO] Found no HLA genes to do LOH analysis")
    print("[INFO] Take a rest now")
    quit(status = 0)
  }

  hla_alleles
}

check_if_alleles_in_ref <- function(alleles, hla_ref_fasta) {
  print(paste(
    "[INFO] Checking if ", paste(alleles, collapse=", "),
    " alleles are defined in the reference",
    sep = ""
  ))

  missing_ones <- alleles[!which(alleles %in% names(hla_ref_fasta))]
  if (length(missing_ones) > 0) {
    print(paste(
      "[WARNING] Found no reference sequence for: ",
      paste(missing_ones, collapse = ","),
      sep = ""
    ))
    print("[WARNING] Trying to find alternatives")
    to_replace <- c()
    for (missing_one in missing_ones) {
      alt <- grep(
        pattern = missing_one, x = names(hla_ref_fasta), value = TRUE
      )[1]
      if (!is.na(alt)) {
        print(paste(
          "[WARNING] Replacing ", missing_one, "with", alt,
          sep = ""
        ))
        allele_type[which(allele_type == missing_one)] <- alt
        to_replace <- c(to_replace, alt)
      } else {
        print(paste(
          "[WARNING] Found no alternative for ", missing_one,
          sep = ""
        ))
        to_replace <- c(to_replace, missing_one)
      }
    }
    # if any of the missing ones still in to_replace, it means they
    # found no alternative from above for loop
    if (any(missing_ones %in% to_replace)) {
      still_missing <- missing_ones[which(missing_ones %in% to_replace)]
      print(paste(
        "[INFO] ", still_missing, " cannot find alternative to replace"
      ))
      print(paste("[INFO] ", alleles, " wont be analyzed for LOH", sep = ""))
      alleles <- NULL
    }
  } else {
    print(paste(
      "[INFO] Found reference sequences for ",
      paste(alleles, collapse = ", "),
      sep = ""
    ))
  }

  alleles
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
  out$a1_aln_start <- a1_aln_start
  out$a1_aln_end <- a1_aln_end
  out$a2_aln_start <- a2_aln_start
  out$a2_aln_end <- a2_aln_end

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

get_seqinfo_from_bam_header <- function(bam) {
  bamf <- BamFile(file = bam)
  bam_header <- scanBamHeader(bamf, what = "targets")
  bam_header$targets
}

check_if_alleles_in_seqinfo <- function(alleles, bam) {
  seqinfo <- get_seqinfo_from_bam_header(bam)
  for (allele in alleles) {
    if (!allele %in% names(seqinfo)) {
      print(paste(
        "[ERROR] ",
        allele,
        " is not in the sequence header of the provided BAM ",
        bam,
        sep = ""
      ))
      quit(status = 1)
    }
  }
}

get_indel_length <- function(cigar) {
  tmp <- unlist(strsplit(gsub("([0-9]+)", "~\\1~", cigar), "~"))
  ins <- grep(pattern = "I", x = tmp)
  del <- grep(pattern = "D", x = tmp)
  total <- sum(as.numeric(tmp[(ins - 1)])) + sum(as.numeric(tmp[del - 1]))
  total
}

#' Function to calculate HLA allele coverage
get_allele_coverage <- function(alleles, bam, hlafa, outdir, sid, min_necnt) {
  allele_coverage <- list()
  outdir <- file.path(outdir, paste(sid, "/", sep = ""))
  parse_dir_path(dir = outdir, create = TRUE)
  for (allele in alleles) {
    seqinfo <- get_seqinfo_from_bam_header(bam = bam)
    bamf <- BamFile(file = bam)

    allele_seq_ln <- seqinfo[which(names(seqinfo) == allele)]
    allele_to_scan <- GenomicRanges::GRanges(
      seqnames = allele,
      ranges = IRanges::IRanges(start = 1, end = allele_seq_ln)
    )
    scan_param <- ScanBamParam(
      flag = scanBamFlag(),
      what = c("qname", "flag", "cigar"),
      which = allele_to_scan,
      tag = "NM"
    )
    obam <- file.path(outdir, paste(sid, allele, "filt.bam", sep = "."))
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
      by=seq_len(nrow(aln_dt))
    ]
    if (nrow(aln_dt) == 0) {
      # this prevents the case in which users provided random HLA alleles
      # FIXME: need to check this at the realignment step
      # a bit too late here
      print(paste(
        "[ERROR] Found no alignments to allele ", allele, " in sample ", sid,
        sep = ""
      ))
      print("[ERRROR] Please check the provided HLA typing results")
      quit(status = 1)
    } else {
      cigar <- NULL # this is to avoid "no visible binding for cigar"
      n_ins <- n_del <- n_mm <- n_ecnt <- NULL
      nread_per_frag <- NULL
      aln_dt[, "n_ins" := length(
        grep(pattern = "I", unlist(strsplit(cigar, "")))
      ), by = seq_len(nrow(aln_dt))]
      aln_dt[, "n_del" := length(
        grep(pattern = "D", unlist(strsplit(cigar, "")))
      ), by = seq_len(nrow(aln_dt))]
      aln_dt[, "n_mm" := nm - apply(.SD, 1, get_indel_length), .SDcols = "cigar"]
      aln_dt[, "n_ecnt" := n_mm + n_ins + n_del]
      aln_dt <- aln_dt[n_ecnt <= min_necnt]
      aln_dt[, "nread_per_frag" := .N, by = "qname"]
      aln_dt <- aln_dt[nread_per_frag == 2]
      filter <- S4Vectors::FilterRules(
        list(function(x) x$qname %in% aln_dt$qname)
      )
      filterBam(bamf, obam, filter = filter, param = scan_param)
      cov_dt <- run_mpileup(bam = obam, hlafa = hlafa, outdir = outdir)
    }

    allele_coverage[[allele]] <- cov_dt
  }

  allele_coverage
}

bin_allele_cov_interval_to_dt <- function(start_pos, end_pos, bin_size) {

  bin_breaks <- seq(start_pos, end_pos, by = bin_size)
  # the last bin can be less than 150, so merged with second last bin
  # end_pos + 2 was from original code
  bin_breaks <- c(bin_breaks[-length(bin_breaks)], end_pos + 2)
  istarts <- bin_breaks[-length(bin_breaks)]
  # this makes sure no end position of last bin not repeating as
  # start position in the next bin
  istarts[2:length(istarts)] <- istarts[2:length(istarts)] + 1
  iends <- bin_breaks[2:length(bin_breaks)]
  bin_dt <- data.table(start = istarts, end = iends)
  bin_dt[, "bin" := paste(start, end, sep = "-")]

  bin_dt
}

makeBins <- function(
  allele, start_pos, end_pos, allele_length, bin_size=150) {

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
    indices <- c(indices, allele_length+1)
  }

  bin_dt <- data.table(
    seqnames=allele, start = istarts, end = iends, bin = indices
  )
  bin_dt[, end := ifelse(end > allele_length, allele_length, end)]

  bin_dt
}

collate_tumor_and_normal_cov_per_allele <- function(
  t_allele_cov_dt, n_allele_cov_dt, allele, multfactor
) {

  setkey(t_allele_cov_dt, start, ref)
  setkey(n_allele_cov_dt, start, ref)
  allele_cov_dt <- t_allele_cov_dt[n_allele_cov_dt]
  allele_cov_dt[, seqnames := ifelse(is.na(seqnames), allele, seqnames)]
  allele_cov_dt[, ":="(end = start, t_dp = dp, n_dp = i.dp)] # nolint
  allele_cov_dt[, ":="(
    i.seqnames = NULL, dp = NULL, i.dp = NULL, qname = NULL, i.qname = NULL,
    flag = NULL, i.flag = NULL
  )]
  #allele_cov_dt[, "logR" := log2((t_dp / n_dp) * multfactor)] # nolint

  allele_cov_dt
}

#call_hla_loh <- function(
#    alleles, tbam, nbam, subj_hla_ref_fasta, outdir,
#    purity, ploidy, multfactor, min_dp, min_necnt,
#    tid = "example_tumor", nid = "example_normal", gamma = 1) {
call_hla_loh <- function(
    row, tbam, nbam, hlaref, outdir,
    purity, ploidy, multfactor, min_dp, min_necnt,
    tid = "example_tumor", nid = "example_normal", gamma = 1) {

  a1 <- row["A1"]
  a2 <- row["A2"]
  alleles_str <- paste(c(a1, a2), collapse = " and ")
  print(paste("[INFO] Analyze LOH for ", alleles_str, sep = ""))

  report <- init_loh_report(a1, a2)

  hla_seq <- read.fasta(hlaref)
  a1_seq <- hla_seq[[a1]]
  a2_seq <- hla_seq[[a2]]
  print(paste(
    "[INFO] Align sequences between ", alleles_str, sep = ""))
  mm <- get_mismatches_bw_alleles(a1_seq, a2_seq)
  report$Num_MM <- length(mm$diffSeq1)

  if (length(mm$diffSeq1) == 0) {
    print(paste(
      "[INFO] there is no difference between the two alleles",
      a1, a2,
      sep = " "
    ))
    print("[INFO] no call will be made. Move to the next HLA gene")
    return(data.table(t(report)))
  }

  if (length(mm$diffSeq1) < 5) {
    print("[WARN] HLA alleles are similar (less than 5 mismatch positions)")
    print("[WARN] Keep that in mind when considering results")
  }

  print(paste(
    "[INFO] Make bins for allele ", a1, sep = ""
  ))
  # FIXME: should simply pass mm variable, rework on getting mm
  a1_bin_dt <- makeBins(
    allele = a1,
    start_pos = mm$a1_aln_start,
    end_pos = mm$a1_aln_end,
    allele_length = length(a1_seq)
  )
  print(paste(
    "[INFO] Make bins for allele ", a1, sep = ""
  ))
  a2_bin_dt <- makeBins(
    allele = a2,
    start_pos = mm$a2_aln_start,
    end_pos = mm$a2_aln_end,
    allele_length = length(a2_seq)
  )

  t_a1_cov <- get_allele_coverage_new(allele = a1, bam = tbam)
  t_a2_cov <- get_allele_coverage_new(allele = a2, bam = tbam)
  n_a1_cov <- get_allele_coverage_new(allele = a1, bam = nbam)
  n_a2_cov <- get_allele_coverage_new(allele = a2, bam = nbam)
  n_a1_cov <- n_a1_cov[count > min_dp]
  n_a2_cov <- n_a2_cov[count > min_dp]
  t_a1_cov <- t_a1_cov[pos %in% n_a1_cov$pos]
  t_a2_cov <- t_a2_cov[pos %in% n_a2_cov$pos]
  if (nrow(n_a1_cov) == 0 || nrow(n_a2_cov) == 0) {
    print("[INFO] Found no coverage for either allele in normal")
    print("[INFO] Move to next HLA gene")
    return(data.table(t(report)))
  }
  print(paste(
    "[INFO] Combine coverage data table for allele ", a1, sep = ""
  ))
  # FIXME: have a prep_allelic_cov function
  a1_cov_dt <- combine_tn_cov(t_dt = t_a1_cov, n_dt = n_a1_cov)
  a1_cov_dt <- bin_allele_cov(cov_dt = a1_cov_dt, bin_dt = a1_bin_dt)
  a1_cov_dt <- estimate_logR(cov_dt = a1_cov_dt, multfactor = multfactor)
  names(a1_cov_dt) <- paste("a1_", names(a1_cov_dt), sep = "")
  a2_cov_dt <- combine_tn_cov(t_dt = t_a2_cov, n_dt = n_a2_cov)
  a2_cov_dt <- bin_allele_cov(cov_dt = a2_cov_dt, bin_dt = a2_bin_dt)
  a2_cov_dt <- estimate_logR(cov_dt = a2_cov_dt, multfactor = multfactor)
  names(a2_cov_dt) <- paste("a2_", names(a2_cov_dt), sep = "")
  print(paste(
    "[INFO] Get median logR for ", alleles_str, sep = ""
  ))
  report$HLA_A1_Median_LogR <- median(a1_cov_dt$a1_logR, na.rm = TRUE)
  report$HLA_A2_Median_LogR <- median(a2_cov_dt$a2_logR, na.rm = TRUE)

  bin_logR_dt <- estimate_binned_logR(
    a1_dt = a1_cov_dt, a2_dt = a2_cov_dt, multfactor = multfactor
  )
  bin_logR_dt[, cn_loss_test_bin := apply(
    .SD, 1, test_cn_loss_using_log_odds_ratio
    ),
    .SDcols=c("a1_bin_t_dp", "a1_bin_n_dp", "a2_bin_t_dp", "a2_bin_n_dp")
  ]
  bin_logR_dt[, capture_bias_bin := a1_bin_n_dp / a2_bin_n_dp]
  report$Num_Bins <- nrow(bin_logR_dt)

  mm_est_dt <- prep_mm_cov(mm = mm, a1_dt = a1_cov_dt, a2_dt = a2_cov_dt)
  if (is.null(mm_est_dt)) {
    print("[INFO] No coverage info found for both alleles at mm sites")
    print("[INFO] Cannot estimate copy numbers. Moving to next HLA gene")
    return(data.table(t(report)))
  }
  mm_est_dt <- estimate_baf(mm_dt = mm_est_dt, bin_dt = bin_logR_dt)

  print("[INFO] Estimate copy number at mismatch sites")
  cn_dt <- estimate_cn(mm_dt = mm_est_dt)
  #na_cn_est <- nb_cn_est <- NA
  print(cn_dt[, c("bin", "logR_combined_bin", "a1_bin_cn", "a2_bin_cn")])
  report$HLA_A1_CN <- median(cn_dt$a1_bin_cn, na.rm = TRUE)
  report$HLA_A2_CN <- median(cn_dt$a2_bin_cn, na.rm = TRUE)

  cn_est_conf <- estimate_cn_conf(cn_dt = cn_dt, which = "a1")
  report$HLA_A1_CN_Lower <- cn_est_conf$est_lower
  report$HLA_A1_CN_Upper <- cn_est_conf$est_upper
  cn_est_conf <- estimate_cn_conf(cn_dt = cn_dt, which = "a2")
  report$HLA_A2_CN_Lower <- cn_est_conf$est_lower
  report$HLA_A2_CN_Upper <- cn_est_conf$est_upper
  print(report)
  stop()

  #print(paste(
  #  "[INFO] Get binned tumor and normal dp for ", alleles_str, sep = ""
  #))
  #a1_cov_dt[, a1_bin_t_dp := as.numeric(median(a1_t_dp, na.rm=TRUE)), by = "a1_bin"]
  #a1_cov_dt[, a1_bin_n_dp := as.numeric(median(a1_n_dp, na.rm=TRUE)), by = "a1_bin"]
  #a2_cov_dt[, a2_bin_t_dp := as.numeric(median(a2_t_dp, na.rm=TRUE)), by = "a2_bin"]
  #a2_cov_dt[, a2_bin_n_dp := as.numeric(median(a2_n_dp, na.rm=TRUE)), by = "a2_bin"]
  #print(a2_cov_dt)

  #print("[INFO] Get local logR corrector")
  #a1_cov_dt[, bin_multfactor := a1_bin_n_dp / a1_bin_t_dp]
  #a2_cov_dt[, bin_multfactor := a2_bin_n_dp / a2_bin_t_dp]
  #local_multfactor <- min(
  #  median(unique(a1_cov_dt, by="a1_bin")$bin_multfactor, na.rm=TRUE),
  #  median(unique(a2_cov_dt, by="a2_bin")$bin_multfactor, na.rm=TRUE)
  #)
  #if (local_multfactor > 1) {
  #  local_multfactor <- 1
  #} else {
  #  if (multfactor > local_multfactor) {
  #    local_multfactor <- multfactor
  #  }
  #}
  #a1_cov_dt[, bin_multfactor := NULL]
  #a2_cov_dt[, bin_multfactor := NULL]

  #a1_cov_dt[, a1_logR := log2(a1_t_dp / a1_n_dp * local_multfactor)]
  #a1_cov_dt[, a1_bin_logR := median(a1_logR, na.rm = TRUE), by=a1_bin]
  #a2_cov_dt[, a2_logR := log2(a2_t_dp / a2_n_dp * local_multfactor)]
  #a2_cov_dt[, a2_bin_logR := median(a2_logR, na.rm = TRUE), by=a2_bin]
  #print(multfactor)
  #print(local_multfactor)

  #print(paste(
  #  "[INFO] Get median logR for ", alleles_str, sep = ""
  #))
  #a1_median_logr <- median(a1_cov_dt$a1_logR, na.rm = TRUE)
  #a2_median_logr <- median(a2_cov_dt$a2_logR, na.rm = TRUE)
  #a1_mm_median_logr <- median(
  #  a1_cov_dt[a1_pos %in% mm$diffSeq1]$a1_logR, na.rm = TRUE
  #)
  #a2_mm_median_logr <- median(
  #  a2_cov_dt[a2_pos %in% mm$diffSeq2]$a2_logR, na.rm = TRUE
  #)
  #print(a1_median_logr)
  #print(a2_median_logr)

  #check_if_alleles_in_seqinfo(alleles, tbam)
  #check_if_alleles_in_seqinfo(alleles, nbam)

  #print(paste(
  #  "[INFO] Get coverage for ", alleles_str, " in tumor", sep = ""
  #))
  #t_hla_cov <- get_allele_coverage(
  #  alleles = alleles,
  #  bam = tbam,
  #  hlafa = subj_hla_ref_fasta,
  #  outdir = outdir,
  #  sid = tid,
  #  min_necnt = min_necnt
  #)
  #print(paste(
  #  "[INFO] Get coverage for ", alleles_str, " in normal", sep = ""
  #))
  #n_hla_cov <- get_allele_coverage(
  #  alleles = alleles,
  #  bam = nbam,
  #  hlafa = subj_hla_ref_fasta,
  #  outdir = outdir,
  #  sid = nid,
  #  min_necnt = min_necnt
  #)
  #n_hla_cov[[a1]] <- n_hla_cov[[a1]][dp > min_dp]
  #t_hla_cov[[a1]] <- t_hla_cov[[a1]][start %in% n_hla_cov[[a1]]$start]
  #n_hla_cov[[a2]] <- n_hla_cov[[a2]][dp > min_dp]
  #t_hla_cov[[a2]] <- t_hla_cov[[a2]][start %in% n_hla_cov[[a2]]$start]

  #if (nrow(n_hla_cov[[a1]]) == 0 || nrow(n_hla_cov[[a2]]) == 0) {
  #  print("[INFO] Found no coverage in normal for either allele")
  #  print("[INFO] Move to next HLA gene")
  #  return(
  #    data.table(
  #      "HLA_A1" = a1, "HLA_A2" = a2,
  #      "HLA_A1_CN" = NA,
  #      "HLA_A1_CN_Lower" = NA,
  #      "HLA_A1_CN_Upper" = NA,
  #      "HLA_A2_CN" = NA,
  #      "HLA_A2_CN_Lower" = NA,
  #      "HLA_A2_CN_Upper" = NA,
  #      "HLA_A1_Median_LogR" = NA,
  #      "HLA_A2_Median_LogR" = NA,
  #      "HLA_A1_MM_Median_LogR" = NA,
  #      "HLA_A2_MM_Median_logR" = NA,
  #      "Median_BAF" = NA,
  #      "HLA_MM_Once_LogR_Paired_Pvalue" = NA,
  #      "HLA_MM_Once_LogR_Unpaired_Pvalue" = NA,
  #      "Num_MM" = length(mm$diffSeq1),
  #      "Num_Bins" = NA,
  #      "Num_CN_Loss_Supporting_Bins" = NA
  #    )
  #  )
  #}

  #print(paste(
  #  "[INFO] Make bins for ", alleles_str, sep = ""
  #))
  #a1_bin_dt <- makeBins(
  #  start_pos = mm$a1_aln_start,
  #  end_pos = mm$a1_aln_end,
  #  allele_length = length(a1_seq),
  #  bin_size = 150
  #)
  #a2_bin_dt <- makeBins(
  #  start_pos = mm$a2_aln_start,
  #  end_pos = mm$a2_aln_end,
  #  allele_length = length(a2_seq),
  #  bin_size = 150
  #)

  # let me not worry about the unique table for now
  #print(paste(
  #  "[INFO] Bin coverage data table for ", alleles_str, sep = ""
  #))
  #a1_cov_dt <- collate_tumor_and_normal_cov_per_allele(
  #  t_allele_cov_dt = t_hla_cov[[a1]],
  #  n_allele_cov_dt = n_hla_cov[[a1]],
  #  allele = a1,
  #  multfactor = multfactor
  #)
  #a1_cov_dt <- bin_allele_cov(
  #  allele_cov_dt = a1_cov_dt,
  #  bin_dt = a1_bin_dt,
  #  allele = a1
  #)
  #names(a1_cov_dt) <- paste(
  #  "a1_", names(a1_cov_dt),
  #  sep = ""
  #)
  #a2_cov_dt <- collate_tumor_and_normal_cov_per_allele(
  #  t_allele_cov_dt = t_hla_cov[[a2]],
  #  n_allele_cov_dt = n_hla_cov[[a2]],
  #  allele = a2,
  #  multfactor = multfactor
  #)
  #a2_cov_dt <- bin_allele_cov(
  #  allele_cov_dt = a2_cov_dt,
  #  bin_dt = a2_bin_dt,
  #  allele = a2
  #)
  #names(a2_cov_dt) <- paste(
  #  "a2_", names(a2_cov_dt),
  #  sep = ""
  #)

  #print(paste(
  #  "[INFO] Get binned tumor and normal dp for ", alleles_str, sep = ""
  #))
  #a1_cov_dt[, a1_bin_t_dp := as.numeric(median(a1_t_dp, na.rm=TRUE)), by = "a1_bin"]
  #a1_cov_dt[, a1_bin_n_dp := as.numeric(median(a1_n_dp, na.rm=TRUE)), by = "a1_bin"]
  #a2_cov_dt[, a2_bin_t_dp := as.numeric(median(a2_t_dp, na.rm=TRUE)), by = "a2_bin"]
  #a2_cov_dt[, a2_bin_n_dp := as.numeric(median(a2_n_dp, na.rm=TRUE)), by = "a2_bin"]

  #print("[INFO] Get local logR corrector ")
  #a1_cov_dt[, bin_multfactor := a1_bin_n_dp / a1_bin_t_dp]
  #a2_cov_dt[, bin_multfactor := a2_bin_n_dp / a2_bin_t_dp]
  #local_multfactor <- min(
  #  median(unique(a1_cov_dt, by="a1_bin")$bin_multfactor, na.rm=TRUE),
  #  median(unique(a2_cov_dt, by="a2_bin")$bin_multfactor, na.rm=TRUE)
  #)
  #if (local_multfactor > 1) {
  #  local_multfactor <- 1
  #} else {
  #  if (multfactor > local_multfactor) {
  #    local_multfactor <- multfactor
  #  }
  #}
  #a1_cov_dt[, bin_multfactor := NULL]
  #a2_cov_dt[, bin_multfactor := NULL]

  #a1_cov_dt[, a1_logR := log2(a1_t_dp / a1_n_dp * local_multfactor)]
  #a1_cov_dt[, a1_bin_logR := median(a1_logR, na.rm = TRUE), by=a1_bin]
  #a2_cov_dt[, a2_logR := log2(a2_t_dp / a2_n_dp * local_multfactor)]
  #a2_cov_dt[, a2_bin_logR := median(a2_logR, na.rm = TRUE), by=a2_bin]

  #print(paste(
  #  "[INFO] Get median logR for ", alleles_str, sep = ""
  #))
  #a1_median_logr <- median(a1_cov_dt$a1_logR, na.rm = TRUE)
  #a2_median_logr <- median(a2_cov_dt$a2_logR, na.rm = TRUE)
  #a1_mm_cov_dt <- a1_cov_dt[a1_start %in% mm$diffSeq1]
  #a2_mm_cov_dt <- a2_cov_dt[a2_start %in% mm$diffSeq2]
  #a1_mm_median_logr <- median(a1_mm_cov_dt$a1_logR, na.rm = TRUE)
  #a2_mm_median_logr <- median(a2_mm_cov_dt$a2_logR, na.rm = TRUE)
  #a1_keep_cols <- c(
  #  "a1_seqnames", "a1_bin", "a1_bin_t_dp", "a1_bin_n_dp", "a1_bin_logR"
  #)
  #a2_keep_cols <- c(
  #  "a2_seqnames", "a2_bin", "a2_bin_t_dp", "a2_bin_n_dp", "a2_bin_logR"
  #)
  #bin_logR_dt <- merge(
  #  unique(a1_cov_dt, by="a1_bin")[, ..a1_keep_cols],
  #  unique(a2_cov_dt, by="a2_bin")[, ..a2_keep_cols],
  #  by.x = c("a1_bin"), by.y = c("a2_bin"),
  #  all.x = TRUE, all.y = TRUE
  #)
  #bin_logR_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
  #bin_logR_dt[, logR_combined_bin := log2((a1_bin_t_dp + a2_bin_t_dp) / (a1_bin_n_dp + a2_bin_n_dp) * local_multfactor)]
  #bin_logR_dt[, cn_loss_test_bin := apply(
  #  .SD, 1, test_cn_loss_using_log_odds_ratio
  #  ),
  #  .SDcols=c("a1_bin_t_dp", "a1_bin_n_dp", "a2_bin_t_dp", "a2_bin_n_dp")
  #]
  #bin_logR_dt[, capture_bias_bin := a1_bin_n_dp / a2_bin_n_dp]

  #print("[INFO] Estimate copy number at mismatch sites")
  #mm_est_dt <- data.table(a1_start = mm$diffSeq1, a2_start = mm$diffSeq2)

  #mm_est_dt <- mm_est_dt[
  #  a1_start %in% a1_cov_dt$a1_start & a2_start %in% a2_cov_dt$a2_start
  #]
  #if (nrow(mm_est_dt) == 0) {
  #  print("[INFO] No coverage info found for both alleles at mm sites")
  #  print("[INFO] Cannot estimate copy numbers. Moving to next HLA gene")
  #  return(
  #    data.table(
  #      "HLA_A1" = a1, "HLA_A2" = a2,
  #      "HLA_A1_CN" = NA,
  #      "HLA_A1_CN_Lower" = NA,
  #      "HLA_A1_CN_Upper" = NA,
  #      "HLA_A2_CN" = NA,
  #      "HLA_A2_CN_Lower" = NA,
  #      "HLA_A2_CN_Upper" = NA,
  #      "HLA_A1_Median_LogR" = a1_median_logr,
  #      "HLA_A2_Median_LogR" = a2_median_logr,
  #      "HLA_A1_MM_Median_LogR" = a1_mm_median_logr,
  #      "HLA_A2_MM_Median_logR" = a2_mm_median_logr,
  #      "Median_BAF" = NA,
  #      "HLA_MM_Once_LogR_Paired_Pvalue" = NA,
  #      "HLA_MM_Once_LogR_Unpaired_Pvalue" = NA,
  #      "Num_MM" = length(mm$diffSeq1),
  #      "Num_Bins" = nrow(
  #        bin_logR_dt[bin > 0 & bin < min(length(a1_seq), length(a2_seq))]
  #      ),
  #      "Num_CN_Loss_Supporting_Bins" = NA
  #    )
  #  )
  #}
  #a1_keep_cols <- c(
  #  "a1_seqnames", "a1_bin", "a1_start", "a1_ref", "a1_t_dp", "a1_n_dp", "a1_logR"
  #)
  #a2_keep_cols <- c(
  #  "a2_seqnames", "a2_bin", "a2_start", "a2_ref", "a2_t_dp", "a2_n_dp", "a2_logR"
  #)
  #mm_est_dt <- merge(mm_est_dt, a1_cov_dt[, ..a1_keep_cols], by="a1_start", all.x=TRUE)
  #mm_est_dt <- merge(mm_est_dt, a2_cov_dt[, ..a2_keep_cols], by="a2_start", all.x=TRUE)
  #mm_est_dt[, ":="(bin = a1_bin, a1_bin = NULL, a2_bin = NULL)]
  #keep_cols <- c(
  #  "bin", "a1_bin_t_dp", "a1_bin_n_dp", "a1_bin_logR",
  #  "a2_bin_t_dp", "a2_bin_n_dp", "a2_bin_logR",
  #  "logR_combined_bin", "capture_bias_bin"
  #)
  #mm_est_dt <- merge(mm_est_dt, bin_logR_dt[, ..keep_cols], by = "bin", all.x=TRUE)

  #mm_est_dt[, baf := a1_t_dp / (a1_t_dp + a2_t_dp)] # nolint
  #mm_est_dt[, baf_combined := baf / capture_bias_bin]
  #mm_est_dt[, "nA_combined_bin" :=
  #  (purity - 1 + baf_combined * 2^(logR_combined_bin / gamma) * # nolint
  #    ((1 - purity) * 2 + purity * ploidy)) / purity,
  #]
  #mm_est_dt[, "nB_combined_bin" :=
  #  (purity - 1 - (baf_combined - 1) * 2^(logR_combined_bin / gamma) * # nolint
  #    ((1 - purity) * 2 + purity * ploidy)) / purity,
  #]
  #mm_est_dt[, ":="(
  #  median_nA_combined_bin = median(nA_combined_bin, na.rm = TRUE),
  #  median_nB_combined_bin = median(nB_combined_bin, na.rm = TRUE)
  #),
  #by = "bin"
  #]

  #mm_est_bin_dt <- unique(mm_est_dt, by = "bin")
  #na_cn_est <- nb_cn_est <- NA
  #na_cn_est <- median(mm_est_bin_dt$median_nA_combined_bin, na.rm = TRUE)
  #nb_cn_est <- median(mm_est_bin_dt$median_nB_combined_bin, na.rm = TRUE)

  #na_cn_pvalue <- nb_cn_pvalue <- 99999.0 
  #na_cn_test <- t_test_with_na(
  # mm_est_bin_dt$median_nA_combined_bin
  #)
  #na_cn_est_conf <- na_cn_test$conf.int
  #na_cn_est_lower <- na_cn_est_conf[1]
  #na_cn_est_upper <- na_cn_est_conf[2]
  #nb_cn_test <- t_test_with_na(
  #  mm_est_bin_dt$median_nB_combined_bin
  #)
  #nb_cn_est_conf <- nb_cn_test$conf.int
  #nb_cn_est_lower <- nb_cn_est_conf[1]
  #nb_cn_est_upper <- nb_cn_est_conf[2]


  #hla_gene <- basename(outdir)
  #out_rds <- file.path(outdir, paste(hla_gene, ".data.rds", sep=""))
  #print(paste("[INFO] Dump intermediate data to file: ", out_rds, sep=""))
  #data_to_save <- list(
  #  mm = mm,
  #  a1_bin_dt = a1_bin_dt,
  #  a2_bin_dt = a2_bin_dt,
  #  a1_cov_dt = a1_cov_dt,
  #  a2_cov_dt = a2_cov_dt,
  #  bin_logR_dt = bin_logR_dt,
  #  mm_est_dt = mm_est_dt
  #)
  #saveRDS(data_to_save, out_rds)

  #names(a1_cov_dt) <- gsub("a1_", "", names(a1_cov_dt))
  #names(a2_cov_dt) <- gsub("a2_", "", names(a2_cov_dt))
  #a1_cov_dt <- normalize_hla_seqname(dt = a1_cov_dt) 
  #a2_cov_dt <- normalize_hla_seqname(dt = a2_cov_dt) 
  #a1_cov_dt[, dp := t_dp]
  #a2_cov_dt[, dp := t_dp]

  #data.table(
  #  "HLA_A1" = a1, "HLA_A2" = a2,
  #  "HLA_A1_CN" = na_cn_est,
  #  "HLA_A1_CN_Lower" = round(na_cn_est_lower, digits = 4),
  #  "HLA_A1_CN_Upper" = round(na_cn_est_upper, digits = 4),
  #  "HLA_A2_CN" = nb_cn_est,
  #  "HLA_A2_CN_Lower" = round(nb_cn_est_lower, digits = 4),
  #  "HLA_A2_CN_Upper" = round(nb_cn_est_upper, digits = 4),
  #  "HLA_A1_Median_LogR" = a1_median_logr,
  #  "HLA_A2_Median_LogR" = a2_median_logr,
  #  "HLA_A1_MM_Median_LogR" = a1_mm_median_logr,
  #  "HLA_A2_MM_Median_logR" = a2_mm_median_logr,
  #  "Median_BAF" = median(mm_est_dt$baf_combined, na.rm = TRUE),
  #  "HLA_MM_Once_LogR_Paired_Pvalue" = paired_test_pval,
  #  "HLA_MM_Once_LogR_Unpaired_Pvalue" = unpaired_test_pval,
  #  "Num_MM" = length(mm$diffSeq1),
  #  "Num_Bins" = nrow(
  #    bin_logR_dt[bin > 0 & bin < min(length(a1_seq), length(a2_seq))]
  #  ),
  #  "Num_CN_Loss_Supporting_Bins" = nrow(
  #    bin_logR_dt[cn_loss_test_bin <= 0.01]
  #  )
  #)
}

count_n_reads_from_bam <- function(
  bam, count_dup = FALSE, only_proper = TRUE, only_primary = TRUE,
  tagfilter = list()
) {
  count_n_reads_df <- countBam(
    bam,
    param = ScanBamParam(flag = scanBamFlag(
      isDuplicate = count_dup,
      isProperPair = only_proper,
      isSecondaryAlignment = !only_primary
    ),
    tagFilter=tagfilter
  ))
  count_n_reads_df$records
}

get_alleles_from_bam <- function(bam) {
  bf <- BamFile(bam)
  # a named integer vector 
  targets <- scanBamHeader(bf)$targets

  sort(names(targets))
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
    by=seq_len(nrow(aln_dt))
  ]
  if (nrow(aln_dt) == 0) {
    print(paste(
      "[ERROR] Found no alignments in the give BAM: ", bam, sep = ""
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
  #aln_dt[, "n_ins" := length(
  #  grep(pattern = "I", unlist(strsplit(cigar, "")))
  #), by = seq_len(nrow(aln_dt))]
  #aln_dt[, "n_del" := length(
  #  grep(pattern = "D", unlist(strsplit(cigar, "")))
  #), by = seq_len(nrow(aln_dt))]
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
    "HLA_A1_CN" = NA,
    "HLA_A1_CN_Lower" = NA,
    "HLA_A1_CN_Upper" = NA,
    "HLA_A2_CN" = NA,
    "HLA_A2_CN_Lower" = NA,
    "HLA_A2_CN_Upper" = NA,
    "HLA_A1_Median_LogR" = NA,
    "HLA_A2_Median_LogR" = NA,
    "HLA_A1_MM_Median_LogR" = NA,
    "HLA_A2_MM_Median_logR" = NA,
    "Median_BAF" = NA,
    "Num_MM" = 0,
    "Num_Bins" = NA,
    "Num_CN_Loss_Supporting_Bins" = NA
  )
}

get_allele_coverage_new <- function(allele, bam) {
  print(paste(
    "[INFO] Get coverage for ", allele, " from ", bam, sep = ""
  ))
  seqinfo <- get_seqinfo_from_bam_header(bam = bam)

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
  pileup_param = PileupParam(
    min_mapq = 20, distinguish_strands = FALSE
  )
  p_dt <- setDT(
    pileup(
      file = bam, scanBamParam = scan_param, pileupParam = pileup_param
  ))
  p_dt[, which_label := NULL]
  p_dt
}

combine_tn_cov <- function(t_dt, n_dt, multfactor) {

  a_t <- unique(t_dt$seqnames)
  a_n <- unique(n_dt$seqnames)
  if ( length(a_t) != 1 || length(a_n) != 1 ) {
    print("[ERROR] Only one seqname is expected in tumor and normal tables")
    quit(status = 1)
  }
  if ( a_t != a_n ) {
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

makeBins <- function(
  allele, start_pos, end_pos, allele_length, bin_size=150) {

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
    indices <- c(indices, allele_length+1)
  }

  bin_dt <- data.table(
    seqnames=allele, start = istarts, end = iends, bin = indices
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
    bin_t_dp = as.numeric(median(t_dp, na.rm=TRUE)),
    bin_n_dp = as.numeric(median(n_dp, na.rm=TRUE)),
    bin_multfactor = median(n_dp / t_dp, na.rm=TRUE)
  ), by="bin"]

  cov_dt
}

estimate_logR <- function(cov_dt, multfactor) {
  cov_dt[, logR := log2(t_dp / n_dp * multfactor)]
  # FIXME: i can also use bin_t_dp and bin_n_dp calculate
  # two should be similar
  cov_dt[, bin_logR := median(logR, na.rm = TRUE), by = "bin"]
  cov_dt
}

estimate_binned_logR <- function(a1_dt, a2_dt, multfactor){

  bin_logR_dt <- merge(
    unique(a1_dt, by="a1_bin"),
    unique(a2_dt, by="a2_bin"),
    by.x = c("a1_bin"), by.y = c("a2_bin"),
    all.x = TRUE, all.y = TRUE
  )
  binned_names <- names(bin_logR_dt)[which(grepl("bin", names(bin_logR_dt)))]
  bin_logR_dt <- bin_logR_dt[, ..binned_names]
  bin_logR_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
  bin_logR_dt[, logR_combined_bin := log2(
    (a1_bin_t_dp + a2_bin_t_dp) /
    (a1_bin_n_dp + a2_bin_n_dp) * multfactor
  )]
  bin_logR_dt

}

prep_mm_cov <- function(mm, a1_dt, a2_dt) {
  mm_est_dt <- data.table(a1_pos = mm$diffSeq1, a2_pos = mm$diffSeq2)
  mm_est_dt <- mm_est_dt[
    a1_pos %in% a1_dt$a1_pos & a2_pos %in% a2_dt$a2_pos
  ]
  if (nrow(mm_est_dt) == 0) { return(NULL) }
  simply_cols <- names(a1_dt)[which(! grepl("bin", names(a1_dt)))]
  simply_cols <- c(simply_cols, "a1_bin")
  mm_est_dt <- merge(
    mm_est_dt, a1_dt[, ..simply_cols], by="a1_pos", all.x=TRUE
  )
  mm_est_dt[, ":="(bin = a1_bin, a1_bin = NULL)]
  if (nrow(mm_est_dt) == 0) { return(NULL) }
  simply_cols <- names(a2_dt)[which(! grepl("bin", names(a2_dt)))]
  mm_est_dt <- merge(
    mm_est_dt, a2_dt[, ..simply_cols], by="a2_pos", all.x=TRUE
  )
  if (nrow(mm_est_dt) == 0) { return(NULL) }

  mm_est_dt
}

estimate_baf <- function(mm_dt, bin_dt) {
  mm_dt[, baf := a1_t_dp / (a1_t_dp + a2_t_dp)] # nolint
  if (! "capture_bias_bin" %in% names(bin_dt)) {
    print("[WARN] Found no column named capture_bias_bin in bin_dt")
    print("[WARN] No capture bias will be corrected for BAF")
    mm_dt[, baf_correct := baf]
    return(mm_dt)
  }
  req_cols <- c("bin", "logR_combined_bin")
  miss_cols <- req_cols[which(! req_cols %in% names(bin_dt))] 
  if ( length(miss_cols) > 0 ) {
    print(paste(
      "[ERROR] Miss ",
      paste(miss_cols, collapse=","),
      " columns in bin_dt", sep = ""
    ))
    print("[ERROR] Cannot continue estimate BAF")
    quit(status= 1)
  }
  mm_dt <- merge(
    mm_dt,
    bin_dt[, c("bin", "logR_combined_bin", "capture_bias_bin")],
    by = "bin",
    all.x=TRUE
  )
  mm_dt[, baf_correct := baf / capture_bias_bin]
  mm_dt

}

estimate_cn <- function(mm_dt) {
  mm_dt[, "a1_cn" :=
    (purity - 1 + baf_correct * 2^(logR_combined_bin / gamma) * # nolint
      ((1 - purity) * 2 + purity * ploidy)) / purity,
  ]
  mm_dt[, "a2_cn" :=
    (purity - 1 - (baf_correct - 1) * 2^(logR_combined_bin / gamma) * # nolint
      ((1 - purity) * 2 + purity * ploidy)) / purity,
  ]
  mm_dt[, ":="(
    a1_bin_cn = median(a1_cn, na.rm = TRUE),
    a2_bin_cn = median(a2_cn, na.rm = TRUE)
  ),
  by = "bin"
  ]

  bin_cn_dt <- unique(mm_dt, by = "bin")

}

estimate_cn_conf <- function(cn_dt, which) {
  col <- NULL
  print(which)
  print(names(cn_dt))
  pattern <- paste(which, "bin_cn", sep="_")
  print(pattern)
  col <- names(cn_dt)[which(grepl(pattern, names(cn_dt)))] 
  if ( is.null(col) || length(col) == 0 ) {
    print(paste(
      "[ERROR] Failed to find cn column in cn_dt for ",
      which, " allele", sep=""
    ))
    quit(status = 1)
  }
  print(col)
  cn_test <- t_test_with_na(cn_dt[[col]])
  cn_est_conf <- cn_test$conf.int
  cn_est_lower <- cn_est_conf[1]
  cn_est_upper <- cn_est_conf[2]
  list(est_lower=cn_est_lower, est_upper=cn_est_upper)

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
  t_rg <- get_read_group_from_bam(bam = args$tbam)
  tid <- get_sm_from_rg(rg = t_rg)
  n_rg <- get_read_group_from_bam(bam = args$nbam)
  nid <- get_sm_from_rg(rg = n_rg)
}

if (args$example) {
  # when using the test data provided by LOHHLA.R set these two lines below
  tid <- "example_tumor"
  nid <- "example_normal"
}

print("[INFO] Getting estimates of ploidy and purity")
tstates_dt <- fread(args$tstates, drop = c(1))
#if (!all(c("ploidy", "purity") %in% names(tstates_dt))) {
if (!all(c("TumorPloidy", "TumorPurityNGS") %in% names(tstates_dt))) {
  print(paste(
    "[ERROR] Miss either purity, or ploidy, or both ",
    "in the CNV model results",
    sep = ""
  ))
  quit(status = 1)
}

# FIXME: when no data provided, assuming ploidy=2, purity=?
ploidy <- tstates_dt[["TumorPloidy"]]
purity <- tstates_dt[["TumorPurityNGS"]]
print("[INFO] Getting estimates of ploidy and purity [DONE]")
print(paste("[INFO] Purity = ", purity, " Ploidy = ", ploidy, sep = ""))
if (is.na(ploidy) || is.na(purity)) {
  print("[INFO] Ploidy and purity estimate are missing from CNV result")
  print("[INFO] Terminate due to missing value to estimate HLA copy number")
  quit(status = 0)
}
print(ploidy)
print(purity)

print("[INFO] Counting sequencing depth from realigned normal BAM")
# FIXME: tagfilter should be generated with args$min_necnt, rather than
# hard-coded
n_seq_depth <- count_n_reads_from_bam(bam = args$nbam, tagfilter=list(NM=c(0,1)))
print("[INFO] Counting sequencing depth from realigned tumor BAM")
t_seq_depth <- count_n_reads_from_bam(bam = args$tbam, tagfilter=list(NM=c(0,1)))
if (t_seq_depth <= 0) {
  print("[ERROR] Found no alignments in the provided tumor BAM")
  quit(status = 1)
}
print(n_seq_depth)
print(t_seq_depth)
multfactor <- n_seq_depth / t_seq_depth
print(paste("[INFO] Normal sequencing depth = ", n_seq_depth, sep=""))
print(paste("[INFO] Tumor sequencing depth = ", t_seq_depth, sep=""))
print(paste("[INFO] Sequencing depth correcting factor = ", multfactor, sep=""))

alleles_n <- get_alleles_from_bam(bam=args$nbam)
alleles_t <- get_alleles_from_bam(bam=args$tbam)
if ( ! all.equal(alleles_n, alleles_t)) {
  print("[ERROR] Different set of alleles detected in normal and tumor BAMs")
  print("[ERROR] Alleles in normal BAM: ", paste(alleles_n, collapse="|"))
  print("[ERROR] Alleles in tumor BAM: ", paste(alleles_t, collapse="|"))
  quit(status = 1)
}

alleles_dt <- data.table(Alleles = alleles_n)
alleles_dt[, HLAGene := tstrsplit(Alleles, "_", keep=2)]
alleles_dt[, HLAGene := paste("hla_", HLAGene, sep="")]
alleles_dt[, Alleles := paste(Alleles, collapse=","), by=HLAGene]
alleles_dt[, c("A1", "A2") := tstrsplit(Alleles, ","), by=HLAGene]
alleles_dt[, Alleles := NULL]
alleles_dt <- unique(alleles_dt, by="HLAGene")

# filter input bams by ecnt
filt_nbam = file.path(args$outdir, paste(nid, ".filt.bam", sep=""))
filt_tbam = file.path(args$outdir, paste(tid, ".filt.bam", sep=""))
if ( ! file.exists(filt_nbam) || ! file.exists(filt_tbam)) {
  filter_bam_by_ecnt(bam = args$nbam, obam=filt_nbam, min_necnt = args$min_necnt)
  filter_bam_by_ecnt(bam = args$tbam, obam=filt_tbam, min_necnt = args$min_necnt)
}

alleles_dt[, apply(
  .SD, 1, call_hla_loh,
  tbam=filt_tbam, nbam=filt_nbam, hlaref=args$hlaref,
  outdir=args$outdir, purity=purity, ploidy=ploidy,
  multfactor=multfactor, min_dp=args$min_cov, min_necnt=args$min_nm,
  tid=tid, nid=nid
), by="HLAGene"]

#hla_genes <- unique(
#  unlist(
#    lapply(
#      strsplit(alleles_n, "_"), function (x) paste(x[1:2], collapse="_")
#)))
#print(hla_genes)


#print("[INFO] Loading HLA allele reference sequences")
#hla_ref_fasta <- read.fasta(args$hlaref)
#print("[INFO] Loading HLA allele reference sequences [DONE]")
#
#hla_alleles_to_analyze <- get_hla_alleles(
#  hla_res = args$hlares, hla_ref_fasta = hla_ref_fasta
#)


#hlaloh_res_list <- list()
#for (i in seq_len(length(hla_alleles_to_analyze))) {
#  hla_gene <- tolower(names(hla_alleles_to_analyze)[i])
#
#  hla_alleles <- hla_alleles_to_analyze[[i]]
#
#  outdir <- file.path(args$outdir, paste(hla_gene, "/", sep = ""))
#  parse_dir_path(dir = outdir, create = TRUE)
#  res_dt <- call_hla_loh(
#    alleles = hla_alleles, tbam = realign_tbam, nbam = realign_nbam,
#    subj_hla_ref_fasta = subj_hla_ref_fasta, outdir = outdir,
#    purity = purity, ploidy = ploidy, multfactor = multfactor,
#    min_dp = args$min_cov, min_necnt = args$min_nm,
#    tid = tid, nid = nid, gamma = gamma
#  )
#  hlaloh_res_list <- append(hlaloh_res_list, list(res_dt))
#}
#
#hlaloh_res_dt <- rbindlist(hlaloh_res_list)
#print(hlaloh_res_dt)
#
#out_res <- file.path(args$outdir, paste(args$subject, ".hlaloh.res.tsv", sep=""))
#fwrite(hlaloh_res_dt, out_res, sep="\t", row.names=FALSE, quote=FALSE)
