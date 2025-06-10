# install.packages('dHSIC')
library(dHSIC)
# install.packages('cdcsis')
library(cdcsis)
# install.packages('foreach')
library(foreach)
# install.packages('doParallel')
library(doParallel)
# functions for fitting models
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("RESIT_fitting.R")

#################################
### Algorithm 1 (MERIT phase 2)
#################################

### Useful functions ###
fit_and_test_independence <- function(x,y,z,alpha,model,parsModel = list(),
                                      indtest, parsIndtest, includeY = 0, 
                                      kernel = F, transform = F){
  # fits x using y and tests against z
  y <- as.matrix(y)
  z <- as.matrix(z)
  mod_fit <- train_model(model,y,x,parsModel,includeY)
  r2 <- mod_fit$residuals
  
  # transform r2 if needed
  if (transform) {
    r2 <- sign(r2 - median(r2))*(x-median(r2))^4
  }
  
  # r2 is numeric, z can be discrete
  Tquan <- indtestAll(fct=indtest,x=z,y=r2,alpha=alpha,pars=parsIndtest,kernel=kernel)
  return(Tquan)
}

indtestAll <- function(fct,x,y,alpha,pars = list(),kernel=F,transform = F)
{
  if (transform){
    y <- sign(y - median(y))*(x-median(y))^4
  }
  result <- fct(X=x,Y=y,alpha=alpha,pars=pars,kernel=kernel)
  return(result)
}

indtestHsic <- function(X,Y,alpha=alpha,pars,kernel=F)
{
  if (kernel){
    return(dhsic.test(X,Y,kernel = c("discrete", "gaussian"), B = 2000, method = "bootstrap"))
  } else{
    return(dhsic.test(Y,X))
  }
}

# distance correlation
indtestdcor <- function(X,Y,alpha=alpha,pars,kernel=F)
{
  return(dcor.test(X,Y, R =1000))
}




### algorithm 1 ###
MERIT <- function(M, alpha = 0.05, alpha_pruning = 0.01, model = train_gam, parsModel = list(), 
                  indtest = indtestHsic, parsIndtest = list(), 
                  force_answer = TRUE, output = TRUE, transform = F, cats_num = 1){
  #M contains the data (each col one component)
  #confounder_check indicates subsets of which size the method tries to omit if it doesn't find any possible sink node
  
  # useful functions
  combine_results <- function(acc, res) {
    if (is.null(acc)) {
      acc <- list(indtest_at_end = res$indtest_at_end, C = res$C, parlen = res$parlen)
    } else {
      acc$indtest_at_end <- c(acc$indtest_at_end, res$indtest_at_end)
      acc$C <- cbind(acc$C, res$C)
      acc$parlen <- c(acc$parlen, res$parlen)
    }
    acc
  }
  
  ### VARIABLE INITIALIZATION ###
  stopping <- 1
  p <- dim(M)[2]-cats_num # the number of variables, excluding Y
  for (i in (p+1) : (p+cats_num)){
    M[,i] <- as.factor(M[,i])
  }
  C <- matrix(0,p+cats_num,p) # for parents, p+1 potential parents * p variables
  err <- matrix(0,p,1) # can be omitted
  S <- 1:p # the set of variables
  par <- matrix(0,p,p + cats_num - 1) # p continuous variables * p (including Y) potential parents
  parlen <- rep(0,p) # number of the parents for the selected variables
  variable <- rep(0,p) # the selected variables
  indtest_at_end <- rep(-1,p) # final independence test CHECK WHETHER TO INCLUDE THE CATS
  includeY <- 1:cats_num # includeY contains the "idx" of the cat variables, ranging from 1 to cats_num
  idYs <- c() # The location of the clusters  
  cat.list <- list() # The cluster of cats removed interactively
  d <- 0
  
  ### CLUSTER INITIALIZATION ###
  cl <- makeCluster(5)
  registerDoParallel(cl)
  
  while(length(S)>1){
    
    # update the pointer
    d <- d+1
    

    
    # compute the pvalues of residuals (with and/or without Y) in parallel
    check <-   foreach(k = 1:length(S), .combine = 'rbind', .packages = c('dHSIC', 'energy'), 
                     .export = c('fit_and_test_independence', 'indtestAll', 
                                 'indtestHsic', 'indtestdcor'), 
                     .inorder = T, .verbose = F, .errorhandling="pass")%dopar%{
                       
                         # Your parallel computation here
                         library(cdcsis)
                         source("RESIT_fitting.R", local = T)
                         print(S)
                         
                         
                         i <- S[k] # the kth component in the remaining set
                         S_new <- S
                         S_new <- S_new[-c(k)] # the predictors
                         
                         Fc.Y <- fit_and_test_independence(M[,i],M[,c(S_new,p+includeY)],M[,c(S_new,p+includeY)],
                                                           alpha,model,parsModel,indtest, parsIndtest, 
                                                           includeY = length(includeY), kernel = F)
                         pvalue <- c(-Fc.Y$p.value) # the negative p-value for the sound model
                         print(pvalue)
                         
                         # print details
                         if(output){
                           if(length(includeY) > 0  & pvalue[1]>-alpha){
                             #    print(paste("Independence rejected: test statistic - critical value =",check[k]))
                             print(paste("Independence with Y rejected: p-value =",-pvalue[1], " variable ", i))
                           }
                           else if (length(includeY) > 0){
                             #    print(paste("Independence not rejected: test statistic - critical value =",check[k]))
                             print(paste("Independence with Y not rejected: p-value =",-pvalue[1], " variable ", i))
                           } 
                         }
                         
                       
                       
                       # the explicit return
                       return(pvalue)
                     } 
    
    print(includeY)
    print(min(check))
    if (length(includeY) > 0 & min(check)>-alpha){
      # this means SOME Y is the only sink at this stage, so we should remove it
      allsubsets <- unlist(lapply((length(includeY) - 1):1, combn, x = includeY, simplify = FALSE), recursive = FALSE)
      if (length(includeY) == 1){
        allsubsets <- allsubsets[-2]
      }
      for (subset.ex in allsubsets){
        check <- foreach(k = 1:length(S), .combine = 'rbind', .packages = c('dHSIC', 'energy'), 
                         .export = c('fit_and_test_independence', 'indtestAll', 
                                     'indtestHsic', 'indtestdcor'), 
                         .inorder = T, .verbose = F)%dopar%{
                           
                           tryCatch({
                             # Your parallel computation here
                             library(cdcsis)
                             source("RESIT_fitting.R", local = T)
                             print(S)
                             
                             
                             i <- S[k] # the kth component in the remaining set
                             S_new <- S
                             S_new <- S_new[-c(k)] # the predictors
                             
                             Fc.Y <- fit_and_test_independence(M[,i],M[,c(S_new,p+subset.ex)],M[,c(S_new,p+subset.ex)],
                                                               alpha,model,parsModel,indtest, parsIndtest, 
                                                               includeY = length(subset.ex), kernel = F)
                             pvalue <- c(-Fc.Y$p.value) # the negative p-value for the sound model
                             
                           }, error = function(e) {
                             cat("Error in iteration", i, ": ", e$message, "\n")
                           })
                           
                           # the explicit return
                           return(pvalue)
                         } 

        if (min(check)<=-alpha | length(includeY) == 1){
          choosedset <- includeY[! includeY %in% subset.ex]
          print("subset")
          print(subset.ex)
          print("choosedset")
          print(choosedset)
          if (length(choosedset) == 0){
            choosedset <- includeY
          }
          bb <- which.min(check)
          idYs <- c(idYs,d) 
          cat.list <- c(cat.list,list(choosedset))
          print("cat.list check")
          print(length(choosedset))
          print(cat.list)
          includeY <- subset.ex
          break
        }
      }

    } else {
      # this means that there is a sink variable
      bb <- which.min(check) # find the one with the least dependence measure
    }
    
    variable[d] <- S[bb] # assigned the selected variable to variable[d], this only includes the continuous ones
    S <- S[-c(bb)]
    parlen[d] <- length(S) + length(includeY) # number of parents for the selected variable
    par[d,1:length(S)] <- S # parent variables, length(S)<=p
    
    if (length(includeY) > 0){
      par[d, (length(S) + 1) : parlen[d]] <- p + includeY 
    }
    if(output){
      print(paste("Possible sink node found:",variable[d]))
      print(paste("causal order (beginning at sink):",paste(variable,collapse=" ")))
    }
  }
  

  
  if(d<p){
    d <- d + 1 # update d
    variable[d] <- S[1] # add the final source node
    # test whether Y is its parent node
    if (length(includeY) > 0){
      allsubsets <- unlist(lapply(length(includeY):1, combn, x = includeY, simplify = FALSE), recursive = FALSE)
      choosedset <- integer(0) # reinitialize the choosedset
      for (subset.ex in allsubsets){
        # MAKE SURE THE KERNEL PARAMETERS ARE CORRECT
        Fc.source <- fit_and_test_independence(M[,S[1]],M[,p+subset.ex],M[,p+subset.ex],alpha,model,
                                               parsModel,indtest, parsIndtest, 
                                               includeY = length(subset.ex), kernel = T, 
                                               transform = transform)
        print("subset")
        print(subset.ex)
        print("includeY")
        print(includeY)
        print(Fc.source$p.value>=alpha)
        if (Fc.source$p.value>=alpha){
          # add Y to its parent set
          parlen[d] <- length(subset.ex)
          choosedset <- includeY[! includeY %in% subset.ex ]
          par[d, 1:parlen[d]] <- p+subset.ex
          print("Y is a potential parent of variable for the last variable")
          includeY <- subset.ex
          if (length(choosedset) > 0){
            idYs <- c(idYs, d) # d = p
            cat.list <- c(cat.list,list(choosedset))
          } else {
            idYs <- c(idYs, d+1) # d = p
            cat.list <- c(cat.list,list(includeY))
            includeY <- integer(0)
          }
          
          break
        }
      }
      # Remaining cats are ordered to the end
      if (length(includeY) > 0){
        if (length(choosedset) == 0){
          idYs <- c(idYs, p) # d = p
        } else {
          idYs <- c(idYs, p+1) # d = p
        }
        cat.list <- c(cat.list,list(includeY))
      } 
    }


  }
  if(output){
    print(paste("causal order (beginning at sink):",paste(variable,collapse=" ")))
    print(paste("removing unnecessary edges..."))
    if(!force_answer){
      print(paste("and performing final independence tests..."))
    }
  }
  print(par)
  print(parlen)
  print(variable)
  print(idYs)
  print(cat.list)
  rm(S)
  # phase 1 done: with treatment on Y        indtest_at_end[d] <- -1
  
  ### parallelized phase 2 
  
  res.phase2 <- foreach(d = 1:p, .combine = combine_results, .init = NULL, 
                        .packages = c('dHSIC', 'energy', 'randomForest'),
                        .export = c('fit_and_test_independence', 'indtestAll',
                                    'indtestHsic', 'indtestdcor'),
                        .inorder = T, .verbose = F) %dopar% {
    print(paste("The current variable is: ",variable[d]))
                          
    source("RESIT_fitting.R", local = T)                 
                          
    C.p2 <- rep(0,p+cats_num)
    
    if(err[d] != 1){
      if (parlen[d] == 0){
        # compute the indtest, parlen and C no need to update!
        indtest_at_end.p2 <- -1
        return(list(indtest_at_end = indtest_at_end.p2, C = C.p2, parlen = 0))
      }
      S<-par[d,1:parlen[d]] 
      for(i in 1:length(S)){
        S_new<-S
        S_new<-S_new[-c(1)] # check if it is okay to remove the 1st variable
        if(length(S)==1){
          tsx <- M[,variable[d]]
          kernel <- (par[d,parlen[d]] > p)
          if (kernel){
            Fc <- indtestAll(indtest,M[,par[d,parlen[d]]],tsx,alpha_pruning,parsIndtest,kernel=kernel,
                             transform = transform)
          } else{
            Fc <- indtestAll(indtest,tsx,M[,par[d,1:parlen[d]]],alpha_pruning,parsIndtest)
          }
          # If there's only one potential parent and no hidden variables
          # we only need to test whether they are dependent or not.
          print(Fc$p.value)
        }
        else{
          includeY <- S_new[which(S_new > p)]
          
          if (length(includeY) > 0){
            factor.position <- which(S_new > p)
            S_n <- c(S_new[-factor.position],includeY) 
          } else {
            S_n <- S_new
          }
          
          Fc <- fit_and_test_independence(M[,variable[d]],M[,S_n],M[,par[d,1:parlen[d]]],
                                          alpha_pruning,model,parsModel,indtest, parsIndtest, 
                                          includeY = length(includeY))  
          print(Fc$p.value)
        }
        if(Fc$p.value>alpha_pruning){ 
          S <- S_new # so it can be removed
        }
        else{
          if(length(S)>1){ # permute the dataset to check the next variable
            tmp <- S[1]
            S[1:(length(S)-1)] <- S[2:length(S)]
            S[length(S)] <- tmp
          }
        }
      }
      
      
       if(!force_answer){
        if(length(S) == 0){
          # so it can happen that all the potential parents are not the actual
          # parents
          tsx <- M[,variable[d]]
          Fc <- indtestAll(indtest,tsx,M[,par[d,1:parlen[d]]],alpha_pruning,parsIndtest)
        }
        else{
          print(S) # TO CHECK THE CLAIM
          kernel <- (length(S) == 1 & S[1] == p+1)
          includeY.final <- (p+1) %in% S
          Fc <- fit_and_test_independence(M[,variable[d]],M[,S],M[,par[d,1:parlen[d]]],
                                          alpha_pruning,model,parsModel,indtest, parsIndtest, 
                                          kernel = kernel, includeY = includeY.final)
        }
        indtest_at_end.p2 <- sign(0.05-Fc$p.value)
        # the final indtest of the variable d
      }
      else{
        #-1: ind., +1 dep.
        indtest_at_end.p2 <- -1
      }
      print(S)
      parlen.p2 <- length(S) # number of the parents of variable d
      C.p2[S] <- rep(1,length(S)) # indicate that the parents of variable d are these S elements
    }
    else{
      #-1: ind., +1 dep.
      indtest_at_end.p2 <- -1
    }
    
    # return: indtest_at_end, C, parlen
    return(list(indtest_at_end = indtest_at_end.p2, C = C.p2, parlen = parlen.p2))
  }
  
  
  indtest_at_end <- res.phase2$indtest_at_end
  C <- res.phase2$C[,order(variable)]
  parlen <- res.phase2$parlen
  
  stopCluster(cl)
  
  if(max(indtest_at_end)<0){
    if(output){
      if(!force_answer){
        print(paste("all independence tests were fine..."))
      }else{
        print(paste("Because of force_answer, the output is always a graph, 
                    and no final ind. test is done to check the quality of the solution."))
      }
    }
    res <- list()
    res$C <- C
    res$idYs <- idYs
    res$order <- variable[length(variable):1]
    res$cat.list <- cat.list
    return(res)
  }
  else{
    print(paste("No solution. Given the proposed order, the final ind. 
                test failed for the following variables:"))
    show(variable[which(indtest_at_end == 1)])
    return(NULL)
  }
}

### algorithm 2 ###

cat_parents_search <- function(res, samples){
  # preprocessing
  graph <- cbind(res$C,matrix(0,nrow = nrow(res$C),ncol = ncol(samples) - ncol(res$C)))
  for (i in 1:length(res$idYs)){
    # i: the index of the cluster of cats
    if (res$idYs[i] > ncol(res$C)){
      return(graph)
    }
    # permute the order of continuous variables
    index <- ncol(res$C) - res$idYs[i] + 1
    potential.parents <- res$order[1:index] # Here, we only have the numerical variables
    n <- length(potential.parents)
    # p <- ncol(samples)  # number of variables
    
    for (cat in res$cat.list[[i]]){
      print(cat)
      graph[res$order[index],cat + ncol(res$C)] <- 1
      if (n == 1){
        next
      }
      cl <- makeCluster(5)
      registerDoParallel(cl)
      # start with the second 
      parents <- foreach(pp = potential.parents[n-1:1], .combine = 'c', 
                         .packages = c('dHSIC'), .inorder = F, 
                         .verbose = F)%dopar%{
                           library(cdcsis)
                           # select the children of these potential parents, and condition on 
                           children <- which(graph[pp,potential.parents]==1)
                           if (length(children) == 0){
                             t <- dhsic.test(samples[,pp],samples[,cat + ncol(res$C)])$p.value
                           } else {
                             t <- cdcov.test(samples[,pp],as.factor(samples[,cat + ncol(res$C)]),samples[,children])$p.value
                           } 
                           if (t <= 0.01){
                             return(pp)
                           }
                         }
      
      stopCluster(cl)
      
      # update the graph
      if (!is.null(parents)){
        graph[parents,cat + ncol(res$C)] <- 1
      }
    }
    
    
  }
  
  
  print("Estimated full causal graph:")
  print(graph)
  return(graph)
}
