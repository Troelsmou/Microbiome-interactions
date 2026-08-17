# Copyright (C) 2026 Troels Mouritzen
# 
# Author: Troels Mouritzen
# Affiliation: Aarhus University
# Email: twm@mbg.au.dk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Description: This script takes a phenotype file, a SNP-association file
# and a genotyping file, and predicts a trait in the phenotyping file
# based on the full set of SNPs, the GWAS associated SNPs of another trait,
# and a random set of SNPs selected from across the genome equal
# in number to the number of GWAS SNPs.

# Associated publication: [PAPER TITLE], [JOURNAL], [YEAR]
# DOI: [DOI OF PAPER]

pacman::p_load(data.table, tidyverse, rrBLUP)

parse_args <- function(args) {
  arg_list <- list()
  i <- 1
  while(i <= length(args)) {
    if(grepl("^--", args[i])) {
      # Remove leading --
      name <- sub("^--", "", args[i])
      # If next arg is not another parameter switch, make it the value
      val <- TRUE
      if(i+1 <= length(args) && !grepl("^--", args[i+1])) {
        val <- args[i+1]
        i <- i + 1
      }
      arg_list[[name]] <- val
    }
    i <- i + 1
  }
  arg_list
}


handle_args <- function(args) {
  args$n_rep = as.numeric(args$n_rep)
  args$n_fold = as.numeric(args$n_fold)
  args$n_random = as.numeric(args$n_random)
  return(args)
}


print_args <- function(args) {
  cat("Arguments:\n")
  for(name in names(args)) {
    cat(paste0(name, ": ", as.character(args[[name]]), "\n"))
  }
  cat("\n")
}


read_files_wrapper <- function(args) {
  files <- read_files(args)
  files <- remove_missing_lines(files)
  files <- get_right_traits(files, args)
  files$GD <- maf_filter(files$GD, threshold = 0.05)
  files$sig <- remove_missing_SNPs(files$GD, files$sig)
  files$pheno <- format_and_remove_covariates(files$pheno)
  files$pheno <- ensure_full_rank(files$pheno)
  return(files)
}


read_files <- function(args) {
  pheno <- read_csv(args$pheno)
  GD <- fread(args$GD) %>%
    as_tibble()
  sig <- read_csv(args$sigfile)
  return(list(pheno = pheno, GD = GD, sig = sig))
}


remove_missing_lines <- function(files) {
  pheno <- files$pheno
  GD <- files$GD
  lines <- intersect(pheno[[1]], GD[[1]])
  pheno <- pheno[pheno[[1]] %in% lines, ]
  GD <- GD[GD[[1]] %in% lines, ]
  files$pheno = pheno
  files$GD = GD
  return(files)
}


get_right_traits <- function(files, args) {
  pheno <- files$pheno
  sig <- files$sig
  if (!args$prediction_trait %in% pheno[[2]]) {
    stop(paste0("Prediction trait ", args$prediction_trait, " not found in phenotype file"))
  }
  pheno <- pheno[pheno[[2]] == args$prediction_trait, ]
  if (!args$gwas_trait %in% sig[[2]]) {
    warning(paste0("GWAS trait ", args$gwas_trait, " not found in significant SNPs file"))
  }
  sig <- sig[sig[[2]] == args$gwas_trait, ]
  files$pheno <- pheno
  files$sig <- sig
  return(files)
}


maf_filter <- function(GD, threshold = 0.05) {
  maf <- colSums(GD[-1]) / (nrow(GD)*2)
  maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  return(GD[c(TRUE, maf > threshold)])
}


remove_missing_SNPs <- function(GD, sig) {
  sig <- sig[sig[[3]] %in% colnames(GD), ]
  return(sig)
}


format_and_remove_covariates <- function(pheno) {
  factors = c()
  singletons = c()
  cols = if (ncol(pheno) >= 4) 4:ncol(pheno) else integer(0)
  for (col in cols) {
    n_levels = length(unique(pheno[[col]]))
    if (n_levels <= 1) {
      singletons = c(singletons, col)
    } else if (class(pheno[[col]]) %in% c("character", "factor")) {
      pheno[[col]] <- as.factor(pheno[[col]])
    }
  }
  if (length(singletons) > 0) {
    pheno <- pheno[, -singletons]
  }
  return(pheno)
}


remove_singleton_covariates <- function(pheno) {
  singletons = c()
  cols = if (ncol(pheno) >= 4) 4:ncol(pheno) else integer(0)
  for (col in cols) {
    n_levels = length(unique(pheno[[col]]))
    if (n_levels <= 1) {
      singletons = c(singletons, col)
    }
  }
  if (length(singletons) > 0) {
    pheno <- pheno[, -singletons]
  }
  return(pheno)
}


ensure_full_rank <- function(pheno) {
  if (ncol(pheno) < 4) {
    return(pheno)
  }
  covariates <- colnames(pheno)[4:ncol(pheno)]
  throwaway <- c()
  for (i in 1:length(covariates)) {
    covariates_loop = if (length(throwaway) > 0) covariates[1:i][-throwaway] else covariates[1:i]
    formula <- as.formula(str_c("~ ", str_c(covariates_loop, collapse = " + ")))
    x = model.matrix(formula, data = pheno)
    rank = qr(x)$rank
    if (rank < ncol(x)) {
      throwaway <- c(throwaway, i)
    }
  }
  if (length(throwaway) > 0) {
    throwaway = throwaway + 3
    pheno <- pheno[, -throwaway]
  }
  return(pheno)
}


make_CVs <- function(pheno, n_folds) {
  for (i in 1:100) {
    lines = unique(pheno[[1]])
    folds = sample(1:n_folds, length(lines), replace = T)
    pheno_folds = tibble(Name = lines, fold = folds)
    colnames(pheno)[1] <- "Name"
    pheno[-2] %>%
      left_join(pheno_folds, by = "Name") %>%
      relocate(fold, .after = "Name") -> pheno_folds
    check <- check_folds(pheno_folds)
    if (check) {
      return(pheno_folds)
    }
  }
  stop("Could not create CVs")
}


check_folds <- function(pheno_folds) {
  if (ncol(pheno_folds) < 4) {
    return(TRUE)
  }
  factors = sapply(pheno_folds[, 4:ncol(pheno_folds)], is.factor) %>%
    which() %>%
    names()
  levels = sapply(pheno_folds[, factors], function(x) length(unique(x)))
  for (f in unique(pheno_folds$fold)) {
    pheno_folds %>%
      filter(fold != f) -> train
    for (i in 1:length(factors)) {
      n_levels = length(unique(train[[factors[i]]]))
      full_rank = ncol(ensure_full_rank(train)) == ncol(train)
      if (n_levels < levels[i] || !full_rank) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

make_fixed <- function(pheno) {
  formula <- if (ncol(pheno) > 3) "~ 0" else "~ 1"
  if (ncol(pheno) > 3) {
    for (i in 4:ncol(pheno)) {
      formula <- str_c(formula, " + ", colnames(pheno)[i])
    }
  }
  formula <- as.formula(formula)
  X = model.matrix(formula, data = pheno)
  return(X)
}


choose_SNPs <- function(sig) {
  sig[[1]] %>%
    unique() -> qtls
  output_snps <- c()
  for (qtl in qtls) {
    snps <- sig[[3]][sig[[1]] == qtl] %>%
      unique()
    snp = NULL
    if (length(snps) > 0) {
      snp <- sample(snps, size = 1)
    }
    output_snps <- c(output_snps, snp)
  }
  return(output_snps)
}

make_K <- function(GD) {
  K = (as.matrix(GD[-1]) %>%
         scale() %>%
         tcrossprod())
  K = K / mean(diag(K))
  rownames(K) <- colnames(K) <- GD$taxa
  return(K)
}

mypredict <- function(test, model) {
  my_fixed = make_fixed(test)
  betas = model$beta
  adjustment = as.numeric(my_fixed %*% betas)
  test %>%
    mutate(gblup = as.numeric(model$u[Name]),
           adjusted = Score - adjustment) -> test
  return(test)
}


cross_validate_gblup <- function(pheno, K) {
  pheno %>%
    filter(fold == max(unique(fold)) + 1) -> out_tibble # empty tibble
  
  for (outfold in unique(pheno$fold)) {
    pheno %>%
      filter(fold != outfold) -> train
    pheno %>%
      filter(fold == outfold) -> test
    
    X = make_fixed(train)
    y = train$Score
    Z = diag(nrow(K))
    colnames(Z) <- rownames(Z) <- colnames(K)
    Z <- Z[train$Name, ]
    
    model <- mixed.solve(K = K, X = X, y = y, Z = Z)
    names(model$u) <- colnames(K)
    
    out_tibble <- bind_rows(out_tibble, mypredict(test, model))
  }
  return(out_tibble)
}

evaluate_prediction_models <- function(sig, pheno, GD, n_folds, n_random) {
  cat("Evaluating prediction models...\n")
  pheno %>%
    make_CVs(n_folds) -> pheno
  
  out_tibble <- tibble(type = c("Full", rep(c("GWAS", "Random"), n_random)),
                       accuracy = rep(NA, n_random * 2 + 1))
  K_full <- make_K(GD)
  
  full_predictions <- cross_validate_gblup(pheno, K_full)
  out_tibble$accuracy[1] <- cor(full_predictions$adjusted, full_predictions$gblup)
  
  for (i in seq_len(n_random)) {
    snps <- choose_SNPs(sig)
    randoms <- sample(colnames(GD[-1]), length(snps), replace = F)
    K_gwas <- make_K(GD[, c("taxa", snps)])
    K_random <- make_K(GD[, c("taxa", randoms)])
    tryCatch({
      gwas_predictions <- cross_validate_gblup(pheno, K_gwas)
      random_predictions <- cross_validate_gblup(pheno, K_random)
      out_tibble$accuracy[i*2] <- cor(gwas_predictions$adjusted, gwas_predictions$gblup)
      out_tibble$accuracy[i*2 + 1] <- cor(random_predictions$adjusted, random_predictions$gblup)
    }, error = function(e) {
      cat("Error in rep ", i, "\n", e$message, "\n")
      out_tibble$accuracy[i*2] <- NA
      out_tibble$accuracy[i*2 + 1] <- NA
    })
    
  }
  return(out_tibble)
}

evaluate_prediction_wrapper <- function(files, args) {
  sig <- files$sig
  pheno <- files$pheno
  GD <- files$GD
  n_folds <- args$n_fold
  n_random <- args$n_random
  n_reps <- args$n_rep
  gwas_trait <- args$gwas_trait
  prediction_trait <- args$prediction_trait
  
  out_tibble <- tibble(type = character(),
                       accuracy = double(),
                       rep = integer())
  for (rep in 1:n_reps) {
    mytime = Sys.time()
    cat("\nRepetition ", rep, " of ", n_reps, " at ", mytime, "\n")
    new_prediction <- evaluate_prediction_models(sig, pheno, GD, n_folds, n_random) %>%
      mutate(rep = rep)
    out_tibble <- bind_rows(out_tibble, new_prediction)
  }
  
  out_tibble %>%
    mutate(gwas_trait = gwas_trait,
           prediction_trait = prediction_trait) %>%
    return()
}

main <- function() {
  args <- commandArgs(T) %>%
    parse_args() %>%
    handle_args()
  print_args(args)
  
  files <- read_files_wrapper(args)
  
  
  evaluate_prediction_wrapper(files, args) %>%
    write_csv(args$output)
}

if (!interactive()) {
  main()
}
