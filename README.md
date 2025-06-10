# Causal discovery in mixed additive noise models

## Mixed-type data Extension for Regression and Independence Testing (MERIT)

This is the implementation of MERIT from the paper:

Yao, R., Verdonck, T., & Raymaekers, J. (2025), [Causal discovery in mixed additive noise models](https://proceedings.mlr.press/v258/yao25a.html), AISTATS2025 (Oral).

Code partially based on the R package [CompareCausalNetworks](https://cran.r-project.org/web/packages/CompareCausalNetworks/index.html)

Required packages:  
install.packages("randomForest")  
install.packages("SID")  
install.packages('dHSIC')  
install.packages('cdcsis')  
install.packages('foreach')  
install.packages('doParallel')  
install.packages('igraph')  

## Code

RESIT_fitting.R: Utilis functions.  
MERIT.R: Implementation of the Mixed-type data Extension for Regression and Independence Testing (MERIT).  
data_simulationcode.R: Utilis functions for generating synthetic datasets.  
data_generation.R: Generating the synthetic datasets for experiments.  
examples.R: Examples for training and evaluation using data from example_data/.  


