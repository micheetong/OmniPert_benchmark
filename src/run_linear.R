# load packages
library(SingleCellExperiment)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
Sys.setenv("BASILISK_EXTERNAL_CONDA"="/Users/mtong/miniconda3")

# define pa
# change variables as needed
pa <- list(
  dataset_name = "southard_rpe1_crispri",   
  test_train_config_path = "~/Desktop/benchmark/southard_RPE1_CRISPRi/southard_data/set2conditions.json",
  pca_dim = 4,
  ridge_penalty = 0.1,
  seed = 1,
  gene_embedding = "training_data",   # or "identity", "zero", "random", or path
  pert_embedding = "training_data",   # or "identity", "zero", "random", or path
  working_dir = "~/Desktop/benchmark/southard_RPE1_CRISPRi",
  result_id = "sep_18_linear",
  out_dir = "~/Desktop/benchmark/southard_RPE1_CRISPRi/southard_data"
)

set.seed(pa$seed)

out_dir <- file.path(pa$out_dir)

# define functions
solve_y_axb <- function(Y, A = NULL, B = NULL, A_ridge = 0.01, B_ridge = 0.01){
  stopifnot(is.matrix(Y) || is(Y, "Matrix"))
  stopifnot(is.null(A) || is.matrix(A) || is(A, "Matrix"))
  stopifnot(is.null(B) || is.matrix(B) || is(B, "Matrix"))
  
  center <- rowMeans(Y)
  Y <- Y - center
  
  if(! is.null(A) && ! is.null(B)){
    stopifnot(nrow(Y) == nrow(A))
    stopifnot(ncol(Y) == ncol(B))
    tmp <- as.matrix(Matrix::solve(t(A) %*% A + Matrix::Diagonal(ncol(A)) * A_ridge) %*% t(A) %*% Y %*% t(B) %*% Matrix::solve(B %*% t(B) + Matrix::Diagonal(nrow(B)) * B_ridge))
  }else if(is.null(B)){
    fit <- lm.fit(A, Y)
    tmp <- as.matrix(Matrix::solve(t(A) %*% A + Matrix::Diagonal(ncol(A)) * A_ridge) %*% t(A) %*% Y)
  }else if(is.null(A)){
    fit <- lm.fit(t(B), t(Y))
    tmp <- as.matrix(Y %*% t(B) %*% Matrix::solve(B %*% t(B) + Matrix::Diagonal(nrow(B)) * B_ridge))
  }else{
    stop("Either A or B must be non-null")
  }
  tmp[is.na(tmp)] <- 0
  list(K = tmp, center = center)
}

# Load data
# Change as needed
folder <- "~/Desktop/benchmark/southard_RPE1_CRISPRi/southard_data" 
sce <- zellkonverter::readH5AD(file.path(folder, "perturb_processed.h5ad"))
set2condition <- rjson::fromJSON(file = file.path(pa$test_train_config_path))
if(! "ctrl" %in% set2condition$train){
  set2condition$train <- c(set2condition$train, "ctrl")
}

# Only keep valid conditions
sce <- sce[,sce$condition %in% unlist(set2condition)]

# Clean up the colData(sce) a bit
sce$condition <- droplevels(sce$condition)
sce$clean_condition <- stringr::str_remove(sce$condition, "\\+ctrl")
training_df <- tibble(training = names(set2condition), condition = set2condition) %>% unnest(condition)
colData(sce) <- colData(sce) %>%
  as_tibble() %>%
  tidylog::left_join(training_df, by = "condition") %>%
  DataFrame()

gene_names <- rowData(sce)[["gene_name"]]
rownames(sce) <- gene_names

baseline <- MatrixGenerics::rowMeans2(assay(sce, "X")[,sce$condition == "ctrl",drop=FALSE])

# Pseudobulk everything
psce <- glmGamPoi::pseudobulk(sce, group_by = vars(condition, clean_condition, training))
assay(psce, "change") <- assay(psce, "X") - baseline

train_data <- psce[,psce$training == "train"]

# Get embeddings
gene_emb <- if(pa$gene_embedding == "training_data"){
  pca <- irlba::prcomp_irlba(as.matrix(assay(train_data, "X")), n = pa$pca_dim)
  rownames(pca$x) <- rownames(train_data)
  pca$x
}else if(pa$gene_embedding == "identity"){
  tmp <- Matrix::Diagonal(n = nrow(psce))
  rownames(tmp) <- rownames(psce)
  tmp
}else if(pa$gene_embedding == "zero"){
  tmp <- Matrix::Diagonal(n = nrow(psce)) - Matrix::Diagonal(n = nrow(psce))
  rownames(tmp) <- rownames(psce)
  tmp
}else if(pa$gene_embedding == "random"){
  tmp <- matrix(rnorm(nrow(psce) * pa$pca_dim), nrow = nrow(psce), ncol = pa$pca_dim)
  rownames(tmp) <- rownames(psce)
  tmp
}else{
  t(as.matrix(data.table::fread(pa$gene_embedding)))
}

pert_emb <- if(pa$pert_embedding == "training_data"){
  pca <- irlba::prcomp_irlba(as.matrix(assay(train_data, "X")), n = pa$pca_dim)
  rownames(pca$x) <- rownames(train_data)
  t(pca$x)
}else if(pa$pert_embedding == "identity"){
  tmp <- Matrix::Diagonal(n = ncol(psce))
  colnames(tmp) <- psce$clean_condition
  tmp
}else if(pa$pert_embedding == "zero"){
  tmp <- Matrix::Diagonal(n = ncol(psce)) - Matrix::Diagonal(n = ncol(psce))
  colnames(tmp) <- psce$clean_condition
  tmp
}else if(pa$pert_embedding == "random"){
  tmp <- matrix(rnorm(ncol(psce) * pa$pca_dim), nrow = pa$pca_dim, ncol = ncol(psce))
  colnames(tmp) <- psce$clean_condition
  tmp
}else{
  as.matrix(data.table::fread(pa$pert_embedding))
}
if(! "ctrl" %in% colnames(pert_emb)){
  pert_emb <- cbind(pert_emb, ctrl = rep(0, nrow(pert_emb)))
}
pert_matches <- match(colnames(pert_emb), train_data$clean_condition)
gene_matches <- match(rownames(gene_emb), rownames(train_data))
if(sum(! is.na(pert_matches)) <= 1){
  stop("Too few matches between clean_conditions and pert_embedding")
}
if(sum(! is.na(gene_matches)) <= 1){
  stop("Too few matches between gene names and gene_embedding")
}

gene_emb_sub <- gene_emb[! is.na(gene_matches),,drop=FALSE]
pert_emb_training <- pert_emb[,! is.na(pert_matches),drop=FALSE]
Y <-  assay(train_data, "change")[na.omit(gene_matches), na.omit(pert_matches),drop=FALSE]
coefs <- solve_y_axb(Y = Y, A = gene_emb_sub, B = pert_emb_training,
                     A_ridge = pa$ridge_penalty, B_ridge = pa$ridge_penalty)

pert_matches_all <- match(psce$clean_condition, colnames(pert_emb))
pert_emb_all <- pert_emb[,pert_matches_all,drop=FALSE]
colnames(pert_emb_all) <- psce$clean_condition

baseline <- baseline[na.omit(gene_matches)]
pred <- as.matrix(gene_emb_sub %*% coefs$K %*% pert_emb_all + coefs$center + baseline)

rownames(pred) <- rownames(psce)[na.omit(gene_matches)]

# Print some nice summary statistics
bind_cols(tibble(cond = psce$clean_condition, training = psce$training), 
          as_tibble(t(assay(psce, "change"))) %>% rename_with(~ paste0("obs-", .x)),
          as_tibble(t(pred - baseline))  %>% rename_with(~ paste0("pred-", .x)))  %>%
  pivot_longer(starts_with(c("pred-", "obs-")), names_sep = "-", names_to = c(".value", "gene")) %>%
  tidylog::drop_na() %>%
  summarize(r2 = cor(obs, pred), .by = c(cond, training)) %>%
  group_by(training) %>%
  skimr::skim(r2)

# Store output
library(rjson)
library(readr) 
out_dir_2 <- "~/Desktop/benchmark/southard_RPE1_CRISPRi/southard_data/linear"
dir.create(out_dir_2, recursive = TRUE, showWarnings = FALSE)

# Convert predictions to a list
pred_list <- as.list(as.data.frame(pred))
write_lines(rjson::toJSON(pred_list), file.path(out_dir_2, "all_predictions.json"))
write_lines(rjson::toJSON(rownames(pred)), file.path(out_dir_2, "gene_names.json"))

#### Session Info
sessionInfo()
