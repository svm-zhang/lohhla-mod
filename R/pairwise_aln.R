
# https://rdrr.io/bioc/Biostrings/src/R/PairwiseAlignments-io.R
# using .pre2postaligned function in writePariwiseAlignments function
extract_aln_in_pos <- function(axset) {
  pos <- seq(Biostrings::start(axset@range), Biostrings::end(axset@range))
  data.table::data.table(
    pos = pos,
    pa_pos = Biostrings:::.pre2postaligned(pos, axset)
  )
}

get_mismatches_bw_alleles <- function(a1, a2, hlaref) {
  hla_seq <- seqinr::read.fasta(hlaref)
  a1_seq <- hla_seq[[a1]]
  a2_seq <- hla_seq[[a2]]
  sigma <- Biostrings::nucleotideSubstitutionMatrix(
    match = 2, mismatch = -1, baseOnly = TRUE
  )
  pair_aln <- Biostrings::pairwiseAlignment(
    Biostrings::DNAString(paste(toupper(a1_seq), collapse = "")),
    Biostrings::DNAString(paste(toupper(a2_seq), collapse = "")),
    substitutionMatrix = sigma, gapOpening = -2, gapExtension = -4,
    scoreOnly = FALSE, type = "local"
  )
  a1_aln <- Biostrings::pattern(pair_aln) # Get the pair_aln for the first sequence
  a2_aln <- Biostrings::subject(pair_aln) # Get the pair_aln for the second sequence

  a1_aln_start <- Biostrings::start(a1_aln)
  a1_aln_end <- Biostrings::end(a1_aln)
  a2_aln_start <- Biostrings::start(a2_aln)
  a2_aln_end <- Biostrings::end(a2_aln)

  p_aln <- extract_aln_in_pos(axset = a1_aln)
  s_aln <- extract_aln_in_pos(axset = a2_aln)
  p_aln <- merge(p_aln, s_aln, by = "pa_pos", all = TRUE)
  p_aln[, ":="(a1_seqnames = a1, a2_seqnames = a2)]
  p_aln[, a1_ref := unlist(strsplit(as.character(a1_aln), split = ""))]
  p_aln[, a2_ref := unlist(strsplit(as.character(a2_aln), split = ""))]
  p_aln <- p_aln[a1_ref != "-" & a2_ref != "-"]
  setnames(p_aln, c("pos.x", "pos.y"), c("a1_pos", "a2_pos"))
  diffSeq1 <- p_aln[a1_ref != a2_ref]$a1_pos
  diffSeq2 <- p_aln[a1_ref != a2_ref]$a2_pos

  if (length(diffSeq1) <= 5) {
    print("[INFO] Found insufficient num mismatches between the two alleles")
    return(NULL)
  }
  list(
    aln_dt = p_aln,
    diffSeq1 = diffSeq1,
    diffSeq2 = diffSeq2,
    a1 = list(start = a1_aln_start, end = a1_aln_end, len = length(a1_seq)),
    a2 = list(start = a2_aln_start, end = a2_aln_end, len = length(a2_seq))
  )
}
