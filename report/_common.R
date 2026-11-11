# Shared setup for every chapter of the book (report/_quarto.yml).
#
# Each chapter renders in its own R session, so the analysis is chained by checkpoints: a chapter
# loads the workspace saved at the end of the chapter before it (objects/checkpoints/NN.RData),
# runs, and saves its own. Render the whole book in order (quarto render report/) after any
# change upstream; a single chapter can be re-rendered on its own once the earlier ones exist.
#
# Interactive use: source("report/_common.R") from the project root, then ckpt_load("NN") to pick up where chapter NN ended.

# Project root = the folder holding .here (Liver4s2_BulkSeq/). Found explicitly because here::here()
# stops at report/, which _quarto.yml marks as a Quarto project.
PROJ_ROOT <- rprojroot::find_root(rprojroot::has_file(".here"))
if (interactive()) setwd(PROJ_ROOT)
knitr::opts_knit$set(root.dir = PROJ_ROOT)             # every path is project-relative
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

suppressPackageStartupMessages({
  library(cqn)        # GC/length normalisation; first, so purrr::map masks mclust::map
  library(tidyverse)
  library(purrr)
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggnewscale)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
})
set.seed(20261105)

CKPT_DIR <- file.path(PROJ_ROOT, "objects", "checkpoints")
ckpt_save <- function(name, env = parent.frame()) {
  dir.create(CKPT_DIR, recursive = TRUE, showWarnings = FALSE)
  save(list = ls(env, all.names = TRUE), envir = env,
       file = file.path(CKPT_DIR, paste0(name, ".RData")))
}
ckpt_load <- function(name, env = parent.frame()) {
  f <- file.path(CKPT_DIR, paste0(name, ".RData"))
  if (!file.exists(f)) stop("Checkpoint ", f, " is missing: render the earlier chapters first.")
  load(f, envir = env)
  invisible(f)
}
