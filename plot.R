#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(argparse)
})

options(width = 600)

grp4_colors <- c("#0072B5FF", "#E18727FF", "#925E9FFF", "#AD002AFF")
a1_color <- grp4_colors[1]
a2_color <- grp4_colors[2]

parse_cmd <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--sample",
    metavar = "STR", type = "character", required = TRUE,
    help = "Specify sample ID/name"
  )
  parser$add_argument("--loh_res",
    metavar = "FILE", type = "character", required = TRUE,
    help = "Specify path to the HLA LOH result"
  )
  parser$add_argument("--loh_dir",
    metavar = "DIR", type = "character", required = TRUE,
    help = "Specify path to the base directory of HLA LOH result"
  )
  parser$parse_args()
}

init_plot_params <- function() {
  plot_theme <- theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "bottom",
    legend.justification = c("center"),
    legend.key = element_blank(),
    legend.text = element_text(size = 12)
  )
  list(theme = plot_theme)
}

plot_cov_profile <- function(
    cov_dt, a1, a2, cov_col, out,
    mask = FALSE, height = 6, width = 12) {
  if (mask) {
    cov_dt[seqnames == a1, seqnames := "Allele1"]
    cov_dt[seqnames == a2, seqnames := "Allele2"]
    a1 <- "Allele1"
    a2 <- "Allele2"
  }
  plot_params <- init_plot_params()
  max_x <- max(cov_dt$pos)
  max_y <- max(cov_dt[[cov_col]])
  m <- ggplot(cov_dt, aes(x = pos, y = cov_dt[[cov_col]])) +
    geom_line(aes(color = seqnames), size = 1.2) +
    scale_color_manual(
      name = "",
      values = c(a1_color, a2_color),
      limits = c(a1, a2)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max_x + 200)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y + 30)) +
    labs(x = "Position", y = "Depth") +
    plot_params$theme
  ggsave(out, m, width = width, height = height)
}

plot_tn_cov_profile <- function(
    cov_dt, a1, a2, out,
    mask = FALSE, height = 6, width = 12) {
  if (mask) {
    cov_dt[seqnames == a1, seqnames := "Allele1"]
    cov_dt[seqnames == a2, seqnames := "Allele2"]
    a1 <- "Allele1"
    a2 <- "Allele2"
  }
  print(nrow(cov_dt))
  plot_params <- init_plot_params()
  max_x <- max(cov_dt$pos)
  max_y <- max(cov_dt$t_dp, cov_dt$n_dp)
  m <- ggplot(cov_dt, aes(x = pos, y = t_dp, color = seqnames)) +
    geom_line(size = 0.8) +
    geom_line(data = cov_dt, aes(x = pos, y = n_dp), linetype = "dotdash") +
    facet_wrap(~seqnames, nrow = 2, scales = "free") +
    scale_color_manual(
      name = "",
      values = c(a1_color, a2_color),
      limits = c(a1, a2)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max_x + 200)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max_y + 30)) +
    labs(x = "Position", y = "Depth") +
    plot_params$theme
  ggsave(out, m, width = width, height = height)
}

plot_logr_profile <- function(
    cov_dt, a1, a2, out,
    mask = FALSE, height = 6, width = 12) {
  if (mask) {
    cov_dt[seqnames == a1, seqnames := "Allele1"]
    cov_dt[seqnames == a2, seqnames := "Allele2"]
    a1 <- "Allele1"
    a2 <- "Allele2"
  }
  a1_median_logr <- median(cov_dt[seqnames == a1]$logR, na.rm = TRUE)
  a2_median_logr <- median(cov_dt[seqnames == a2]$logR, na.rm = TRUE)
  plot_params <- init_plot_params()
  max_x <- max(cov_dt$pos)
  min_y <- min(cov_dt$logR)
  max_y <- max(cov_dt$logR)
  m <- ggplot(cov_dt, aes(x = pos, y = logR)) +
    geom_line(aes(color = seqnames), size = 1.2) +
    geom_hline(
      aes(yintercept = a1_median_logr),
      color = a1_color, linetype = "dashed"
    ) +
    geom_hline(
      aes(yintercept = a2_median_logr),
      color = a2_color, linetype = "dashed"
    ) +
    scale_color_manual(
      name = "",
      values = c(a1_color, a2_color),
      limits = c(a1, a2)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max_x + 200)) +
    scale_y_continuous(expand = c(0, 0), limits = c(min_y - 1, max_y + 1)) +
    labs(x = "Position", y = "LogR") +
    plot_params$theme
  ggsave(out, m, width = width, height = height)
}

plot_baf_profile <- function(mm_dt, out, mask = FALSE, height = 6, width = 12) {
  if (mask) {
    mm_dt[, a1_seqnames := "Allele1"]
  }
  a1 <- unique(mm_dt$a1_seqnames)
  plot_params <- init_plot_params()
  max_x <- max(mm_dt$a1_pos)
  m <- ggplot(data = mm_dt, aes(x = a1_pos, y = baf_correct)) +
    geom_point(aes(color = a1_seqnames), size = 0.5, stroke = NA) +
    geom_hline(aes(yintercept = 0.5), color = "grey", linetype = "dashed") +
    scale_color_manual(
      name = "",
      values = c(a1_color),
      limits = c(a1)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max_x + 50)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
    labs(x = "Position", y = "B-allele Frequency") +
    plot_params$theme
  ggsave(out, m, width = width, height = height)
}

plot_profile <- function(hla_gene, dat, plot_dir) {
  mm <- dat$mm
  cov_dt <- dat$cov_dt
  alleles <- unique(cov_dt$seqnames)
  a1 <- alleles[1]
  a2 <- alleles[2]
  a1_aln_start <- mm$a1$start
  a2_aln_start <- mm$a2$start
  cov_dt[seqnames == a1, pos := pos + a2_aln_start - 1]
  cov_dt[seqnames == a2, pos := pos + a1_aln_start - 1]
  out <- file.path(plot_dir, paste(hla_gene, ".t_dp.pdf", sep=""))
  plot_cov_profile(
   cov_dt = cov_dt, a1 = a1, a2 = a2, cov_col = "t_dp", out = out
  )
  out <- file.path(plot_dir, paste(hla_gene, ".n_dp.pdf", sep=""))
  plot_cov_profile(
   cov_dt = cov_dt, a1 = a1, a2 = a2, cov_col = "n_dp", out = out
  )
  out <- file.path(plot_dir, paste(hla_gene, ".tn_dp.pdf", sep=""))
  plot_tn_cov_profile(
    cov_dt = cov_dt, a1 = a1, a2 = a2, out = out
  )
  out <- file.path(plot_dir, paste(hla_gene, ".logR.pdf", sep=""))
  plot_logr_profile(
    cov_dt = cov_dt, a1 = a1, a2 = a2, out = out
  )
  out <- file.path(plot_dir, paste(hla_gene, ".baf.pdf", sep=""))
  plot_baf_profile(mm_dt = dat$mm_dt, out = out)
}

plot_hlaloh <- function(dt, sample, wkdir) {
  hla_gene <- dt[["HLAGene"]]

  hla_gene_dat_file <- file.path(wkdir, paste(hla_gene, ".rds", sep = ""))
  if (!file.exists(hla_gene_dat_file)) {
    print(paste(
      "[INFO] Found no intermediate data for HLA gene: ",
      hla_gene_dat_file,
      sep = ""
    ))
    print("[INFO] Move next HLA gene to plot")
    return(NULL)
  }
  hla_gene_dat <- readRDS(hla_gene_dat_file)
  contents <- ls(hla_gene_dat)
  if (length(contents) == 0) {
    print(paste(
      "[INFO] Found no content in HLA gene data: ",
      hla_gene_dat_file,
      sep = ""
    ))
    print("[INFO] Move next HLA gene to plot")
  }

  plot_dir <- file.path(wkdir, paste(sample, "_plots", sep = ""))
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  plot_profile(hla_gene, hla_gene_dat, plot_dir)
  return(NULL)
}

run_plot <- function() {
  args <- parse_cmd()

  loh_res_dt <- fread(args$loh_res)

  loh_res_dt <- loh_res_dt[, apply(
    .SD, 1, plot_hlaloh,
    sample = args$sample,
    wkdir = args$loh_dir
  )]
}

run_plot()
