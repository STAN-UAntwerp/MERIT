### Test script fot algorithm 1 and 2 ###
#install.packages("randomForest")
library(randomForest)
#install.packages("SID")
library(SID)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("MERIT.R")
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])


RecallPrecision <- function(g1, g2){
  tp <- length(which((g1 == 1) & (g2 == 1)))
  fp <- length(which((g1 == 0) & (g2 == 1)))
  p <- sum(g1)
  p_e <- sum(g2)
  pr <- list(Precision = tp/p_e, Recall = tp/p)
  return(pr)
}

evaluate_wrapper <- function(g1, g2){
  # g1 and g2 are matrices
  SHD <- hammingDist(g1, g2)
  SID <- structIntervDist(g1, g2)
  pr <- RecallPrecision(g1, g2)
  f1 <- (2*pr$Precision*pr$Recall)/(pr$Precision+pr$Recall)
  print(paste("SHD:", SHD))
  print("SID:")
  print((SID))
  print(pr)
  print("F1 score:")
  print(f1)
  result_df <- data.frame("SHD" = SHD, "SID"=SID, "F1"=f1)
}
dir.create(file.path("estimated/"), showWarnings = FALSE)

# algorithm wrapper, from data generation to evalutation
algo_wrapper_rep <- function(indtest = indtestHsic, transform = F, 
                             dag_type = 'linear', p, noise.type, g.type, 
                             function.type, rep){
  # data_type: linear, or mlp
  
  # fetch data
  if (dag_type == 'linear'){
    dataset <- readRDS(paste0("example_data/data_linear_logis_noise_example.rds"))
    print(dataset$graph)
    p <- nrow(dataset$graph)
    # compute results
               res <- MERIT(M = as.matrix(dataset$samples), model = train_linear, 
               force_answer = T, indtest = indtest, cats_num = floor(p/2))

  } else {
    dataset <- readRDS(paste0("example_data/data_mlp_uniform_example.rds"))
    print(dataset$graph)
    p <- nrow(dataset$graph)
    res <- MERIT(M = as.matrix(dataset$samples), model = train_randomforest, 
                 force_answer = T, indtest = indtest, cats_num = floor(p/2))
    }
    
  print(res)
  g_estimated <- cat_parents_search(res, dataset$samples)
  g_truth <- dataset$graph
  
  # evaluation
  evaluate_wrapper(g_truth, g_estimated)
  }
  

algo_wrapper_rep(dag_type = 'linear')
algo_wrapper_rep(dag_type = 'mlp')
