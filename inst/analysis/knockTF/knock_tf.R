#!/usr/bin/env R
devtools::load_all()
remap_promoters_f5 <- xcoredata::remap_promoters_f5()
promoters_f5_core <- xcoredata::promoters_f5_core()

# UTILS ========================================================================
topOrSigRes <- function(res) {
  nsig <- sum(abs(res[["z_score"]]) > 1.95)
  if (nsig >= 10) {
    head(res, nsig)
  } else {
    head(res, 10)
  }
}

getRank <- function(res, TF) {
  m <- grep(TF, res$name)
  ifelse(length(m), m[1], 0)
}

baseName <- function(x, suffix = "") {
        x %>% 
        sub(pattern = ".*/", replacement = "") %>%
        sub(pattern = suffix, replacement = "")
}

symbol2fantom <- setNames(promoters_f5_core$name, promoters_f5_core$SYMBOL)

processCounts <- function(in_counts) {
  counts_symbol <- in_counts[, -1] %>% as.matrix() %>% type.convert(as.is = TRUE)
  rownames(counts_symbol) <- in_counts[, 1]
  counts_symbol <-
    counts_symbol[rownames(counts_symbol) %in% promoters_f5_core$SYMBOL, ]
  rownames(counts_symbol) <- symbol2fantom[rownames(counts_symbol)]
  
  return(counts_symbol)
}

makeDesign <- function(mat) {
  cond <-
    ifelse(grepl(pattern = "k", x = colnames(mat)),
           yes = "knockdown",
           no = "control")
  ncond <- cond %>% unique() %>% length()
  design <- diag(ncond)[cond %>% factor() %>% as.numeric(), ]
  rownames(design) <- colnames(mat)
  colnames(design) <- cond %>% factor() %>% levels()
  
  return(design)
}

makeMAE <- function(counts_symbol, ex_design) {
  design_wo_base <- ex_design
  design_wo_base[, "control"] <- 0
  mae <- MultiAssayExperiment(
    experiments = MultiAssayExperiment::ExperimentList(
      U = matrix(
        data = rowMeans(counts_symbol[, ex_design[, "control"] == 1]),
        ncol = 1,
        dimnames = list(rownames(counts_symbol), "u")
      ),
      Y = counts_symbol[, ex_design[, "control"] == 0]
    ),
    metadata = list(design = design_wo_base)
  )
  
  return(mae)
}

# ==============================================================================

path_to_files <- system.file("inst", "extdata", "knockTF", package = "xcore")
dump_isilon <- "/mnt/restricted/isihome/mmigdal/projects/knockTF/"

argv <- commandArgs(trailingOnly=TRUE) %>% as.numeric()

results <- list()
in_files <-
  list.files(
    path = path_to_files,
    pattern = "DataSet_02.*.txt",
    full.names = TRUE
  )
samples <- baseName(in_files, ".txt")

doMC::registerDoMC(cores = 2L) # usually there are just two knockdown replicates
for (i in seq(argv[1], argv[2])) {
  set.seed(314159265)
  in_counts <- read.table(file = in_files[i], header = TRUE)
  counts_symbol <- processCounts(in_counts)
  ex_design <- makeDesign(counts_symbol)
  mae <- makeMAE(counts_symbol, ex_design)
  
  mae <- addSignatures(mae, remap = remap_promoters_f5)
  mae <- filterSignatures(mae, min = 0.05, max = 0.95)
  
  res <- modelGeneExpression(mae = mae,
                             xnames = c("remap"),
                             nfolds = 10)
  save(res, file = paste0(samples[i], ".rda"))
  cmd_mv_to_isilon <- paste0("mv", " ", samples[i], ".rda", " ", dump_isilon)
  system(cmd_mv_to_isilon, wait = FALSE)
  
  results[[samples[i]]] <- res$results$remap
}

save(results, file = paste0("results_", argv[1], "_", argv[2], ".rda"))
