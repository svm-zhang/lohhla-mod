require(Rsamtools)
require(data.table)

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
  if (is.null(rg)) {
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

count_n_reads_from_bam <- function(bam, which, tagfilter = list()) {
  seqinfo <- extract_seqinfo_from_bam(bam = bam)
  allele_seq_ln <- seqinfo[which(names(seqinfo) %in% which)]
  dt <- data.table(seqnames = names(allele_seq_ln), end = allele_seq_ln)
  dt[, start := 1]
  allele_to_scan <- GenomicRanges::makeGRangesFromDataFrame(dt)
  count_n_reads_df <- countBam(
    bam,
    param = ScanBamParam(
      which = allele_to_scan,
      flag = scanBamFlag(
        isDuplicate = FALSE,
        isProperPair = TRUE,
        isNotPassingQualityControls = FALSE,
        isSupplementaryAlignment = FALSE,
        isSecondaryAlignment = FALSE
      ),
      tagFilter = tagfilter
    )
  )
  sum(count_n_reads_df$records)
}

get_indel_length <- function(cigar) {
  tmp <- unlist(strsplit(gsub("([0-9]+)", "~\\1~", cigar), "~"))
  ins <- grep(pattern = "I", x = tmp)
  del <- grep(pattern = "D", x = tmp)
  total <- sum(as.numeric(tmp[(ins - 1)])) + sum(as.numeric(tmp[del - 1]))
  total
}

filter_bam_by_ecnt <- function(bam, obam, min_necnt = 1) {
  bamf <- BamFile(file = bam)

  # including secondary alignments, otherwise, one random
  # allele out of 2 has more coverage than expected, while
  # the other has lower-than-expected coverage
  scanflag <- scanBamFlag(
    isProperPair = TRUE, isSecondaryAlignment = NA, isDuplicate = FALSE,
    isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE
  )
  scan_param <- ScanBamParam(
    flag = scanflag,
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
  # because we are including secondary alignments, here using mod operator
  aln_dt <- aln_dt[nread_per_frag %% 2 == 0]
  filter <- S4Vectors::FilterRules(
    list(function(x) x$qname %in% aln_dt$qname)
  )
  obam <- filterBam(bamf, obam, filter = filter, param = scan_param)
}

extract_allele_coverage <- function(allele, bam, hlaref, min_dp = 0) {
  dt <- fread(text = system2(
    command = "samtools",
    args = c(
      "mpileup", "-f", hlaref,
      "--rf", 2, "--ff", 3584, "-Q", 20, "-r", allele, bam
    ),
    stdout = TRUE,
    stderr = FALSE,
    wait = TRUE
  ), select = c(1, 2, 3, 4))
  # FIXME: should I return NULL?
  if (nrow(dt) == 0) {
    print(paste(
      "[ERROR] No pileup generated for allele ",
      allele,
      " from bam ",
      bam,
      sep = ""
    ))
    quit(status = 1)
  }
  names(dt) <- c("seqnames", "pos", "nucleotide", "count")
  dt <- dt[count > min_dp]
  dt
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
  scanflag <- scanBamFlag(
    isProperPair = TRUE, isSecondaryAlignment = FALSE, isDuplicate = FALSE,
    isNotPassingQualityControls = FALSE, isSupplementaryAlignment = FALSE
  )
  scan_param <- ScanBamParam(
    flag = scanflag,
    mapqFilter = 20,
    which = allele_to_scan,
  )
  pileup_param <- PileupParam(
    min_mapq = 20, distinguish_strands = FALSE, max_depth = 9999
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
