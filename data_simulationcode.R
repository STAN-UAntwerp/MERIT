
library(igraph)

### useful functions ###
# basic functions
runif2 <- function(n){
  sign <- rbinom(n,1,0.5)
  x.p <- runif(sum(sign), 0.5, 2)
  x.n <- runif(n-sum(sign), -2, -0.5)
  x <- rep(0,n)
  x[sign == 0] <- x.n
  x[sign == 1] <- x.p
  return(x)
}

indexy <- function(p){
  x <- rep(0,p)
  id <- sample(1:p, 1)
  x[id] <- 1
  return(x)
}

### simulation function for linear connections ###
compute_samples <- function(coef.m, noise, idxy, res, n_level = 1, cats_num = 1,
                            pre_cat = NULL, coef_cat = NULL, child_cat = NULL,
                            coef.m.c = NULL, intercept_cat = NULL){
  child_cat_all <- unique(unlist(child_cat))
  # generate data given the graph and the noise
  for (j in 2: ncol(coef.m)){
    if (j %in% idxy){
      
      id_j <- which(idxy == j)
      
      if (n_level > 1){
        # print(as.vector(pre_cat + (coef.m.c%*%res)))
        res[j] <- which.max(as.matrix(pre_cat)[id_j,] + (coef.m.c[[id_j]]%*%res)) 
      } else {
        # This is not used
        res[j] <- sign(sum(coef.m[,j]*res) + noise[j])
      }
      
    } else {
      res[j] <- sum(coef.m[,j]*res) + noise[j]
      for (i in 1:cats_num){
        if (j %in% child_cat[[i]]){
          res[j] <- res[j] + intercept_cat[[i]][which(child_cat[[i]]==j),res[idxy[i]]]
        }
      }
      
    } 
  }
  return(res)
}



sim.lingam2 <- function(n, p, ify = T, coef.m = NULL, m = NULL, idxy = NULL, 
                        y_no_parent= F, n_level = 1, g.type = "ER", 
                        noise.type = "uniform", cats_num = 1){
  
  # n number of samples, p number of predictors
  
  # simulate the graph
  if (is.null(m)){
    if (g.type == 'ER'){
      rdb <- rbinom(p*(p-1)/2, 1, 2/(p-1))
      m <- diag(x=0, p)
      upperm = upper.tri(m, diag = F)
      m[upperm] = rdb
    } else if (g.type == 'SF'){
      m.p <- 2
      g <- barabasi.game(n = p, m = m.p, directed = T)
      m <- as.matrix(as_adjacency_matrix(g))
      m <- t(m)
    }
  } 
  
  # simulate the coefficients
  if (is.null(coef.m)){
    coeff <- runif2(p*(p-1)/2)
    coef.m <- diag(x=0, p)
    upperm.coef <-  upper.tri(coef.m, diag = F)
    coef.m[upperm.coef] <-  coeff
  } else {
    # re-permute the matrix if input, similar to m
    if (idxy == 1){
      coef.m <- rbind(coef.m[p,], coef.m[1:(p-1), ])
      coef.m <- cbind(coef.m[,p], coef.m[, 1:(p-1)])
    } else if (idxy != p) {
      coef.m <- rbind(coef.m[1:(idxy-1),], coef.m[p,], coef.m[idxy:(p-1), ])
      coef.m <- cbind(coef.m[,1:(idxy-1)], coef.m[,p], coef.m[, idxy:(p-1)])
    } 
  }
  
  print(coef.m)
  print(m)
  coef.m <- coef.m*m
  
  # generate noises
  if (noise.type == "uniform"){
    noise <- matrix(runif(n*p, -3, 3), p, n) # simple noise 
  } else if (noise.type == "logis"){
    noise <- matrix(rlogis(n*p, 0, 2), p, n)
  }
  
  # sample position for Y
  if (ify){
    idxy <- sample(1:p, cats_num)
    # adjust the noise for Y, let's say NORMAL distribution
    # not used if n_level > 1
    for (i in idxy){
      noise[i,] <- rnorm(n, 0, 1)
    }
  } else if (is.null(idxy)){
    # THIS CORRESPONDS TO THE CASE WHEN THERE'S NO Y
    idxy <- -1
  } else {
    # for reproduction, idxy is specified
    idxy <- idxy
    noise[idxy,] <- rnorm(n, 0, 1)
  }
  print("idxy")
  print(idxy)
  res <- matrix(0, p, n)
  
  
  if (cats_num > 1){
    for (i in idxy){
      for (j in idxy[2:cats_num]){
        m[i, j] <- 0
        m[j, i] <- 0
        coef.m[i, j] <- 0
        coef.m[j, i] <- 0
      }
    }
  }
  
  child_cats <- list()
  parent_cats <- list()
  for (i in 1:cats_num){
    child_cats[[i]] <-which(m[idxy[i],] == 1)
    parent_cats[[i]] <- which(m[,idxy[i]] == 1)
  }

  
  if (n_level>1) {
    coef.m.c.list <- list()
    intercept_cat.list <- list()
    unscaled_prob_cats <- array(0, dim = c(cats_num, n, n_level))
    for (i in 1:cats_num){
      current_id <- idxy[i]
      shift_cat <- sample(seq(-2,2,0.05), n_level, replace = T)
      print("shift_cat")
      print(shift_cat)
      # sample logis noise
      noise_cat <- matrix(rlogis(n*n_level), nrow = n)
      unscaled_prob_cats[i, , ] <- sweep(noise_cat, 2, shift_cat, FUN = "+") 
      if (length(parent_cats[[i]]) == 0){
        coef.m.c <- matrix(0, nrow = n_level, ncol = p)
      } else {
        coef.m.c <- matrix(0, nrow = n_level, ncol = p)
        coef.m.c[,parent_cats[[i]]] <- runif2(length(parent_cats[[i]])*n_level)
      }
      print("coef.m.c")
      print(coef.m.c)
      # sample coefficient from cat
      coef.m[current_id,] <- 0
      intercept_cat <- matrix(runif(length(child_cats[[i]])*n_level, -2, 2), nrow = length(child_cats[[i]]))
      print("Intercept_cat:")
      print(intercept_cat)
      coef.m.c.list[[i]] <- coef.m.c
      intercept_cat.list[[i]] <- intercept_cat
    }
  
  } else {
    unscaled_prob_cat <- NULL
    intercept.cat <- NULL
    coef.m.c <- NULL
  }
  
  # adjust res and causal graph
  for (i in 1:cats_num){
    if (idxy[i] == 1){
      if (n_level == 1){
        res[1,] <- sign(noise[1,])
      } else {
        # sample shifts
        res[1,] <- apply(unscaled_prob_cats[i,,],1,which.max)
      }
      break 
    }
    else{
      res[1,] <- noise[1,]
    }
  }
  
  # compute the models
  # pre_Cat n x level
  print(dim(unscaled_prob_cats[,i,]))
  samples <- do.call(rbind, 
                     lapply(1:n, function(i){
                       compute_samples(coef.m, noise[,i], idxy, res[,i], 
                                       n_level = n_level, cats_num = cats_num, 
                                       pre_cat = unscaled_prob_cats[,i,],
                                       child_cat = child_cats, 
                                       intercept_cat = intercept_cat.list,
                                       coef.m.c = coef.m.c.list)}))
  
  # permutation to make Y the last variable
  print("samples before permutation")
  print(samples)
  
  ordered_idxy <- idxy[order(idxy)]
  for (i in 1:cats_num){
    current_id <- ordered_idxy[i] - (i - 1) 
    print(current_id)
    
    m <- rbind(m[-current_id,], m[current_id,])
    coef.m <- rbind(coef.m[-current_id,], coef.m[current_id, ])
    m <- cbind(m[,-current_id], m[,current_id])
    coef.m <- cbind(coef.m[,-current_id], coef.m[,current_id])  
    
    # also samples and noise!
    noise <- rbind(noise[-current_id, ], noise[current_id, ])
    samples <- cbind(samples[, -current_id], samples[, current_id])
  }
  
  # generate Y
  return(list(graph = m, coef = coef.m, noise = noise, 
              samples = samples, idxy = idxy, coef.m.c = coef.m.c.list))
}


### simulation function for MLP connections ###

# generation of a two-layer MLP
mlp <- function(X, size, seed, shift, onehot = NULL, cat_level = 1){
  
  X_OUTPUT <- NULL
  # TO GENERATE CATEGORICAL DATA, WE NEE TO COMPUTE THE COORDINATES
  for (rep in 1:cat_level) {
    if (!is.null(seed)){
      set.seed(seed[rep]) 
    }
    
    # for source nodes
    if (sum(X) == 0){
      # return a constant value
      X_OUTPUT <- rbind(X_OUTPUT, rep(0,nrow(X)))
      next
    }
    
    # onehot encoding and dimension decision
    if (!is.null(onehot)){
      for (cat_id in onehot){
        cat <- as.factor(X[,cat_id])
        print(length(levels(cat)))
        one_hot_encoded <- do.call(cbind, lapply(1:cat_level, 
                                                 function(i){res <-  rep(0,nrow(X))
                                                 res[X[,onehot]==i] <- 1
                                                 return(res)}))
        X <- cbind(X,one_hot_encoded) 
      }
      X <- X[,-onehot] 
    }
    dim <-  ncol(X) # updated dimension
    
    # simulate matrix
    M_input = matrix(rnorm(dim*size), nrow = size)
    
    X_inter = t(apply(X, 1, function(x){M_input%*%(x)}))*4
    
    X_inter = plogis(X_inter) # n x size
    hist(X_inter)
    M_output = rnorm(size)
    X_output = apply(X_inter, 1, function(x){M_output%*%(x)}) # 1 x n
    
    # adjust the output
    if (is.null(shift)){
      bias_output <- 0
    } else {
      bias_output <- shift
    }
    print(bias_output)
    X_output <- X_output - bias_output
    X_OUTPUT <- rbind(X_OUTPUT, X_output) 

  }
  
  if (cat_level>1) {
    X_OUTPUT <- X_OUTPUT+1
  }
  return(X_OUTPUT) # n_level X n or 1 x n
  
}


compute_samples_mlp <- function(noise, idxy, res, seeds, m, size = 200, 
                                shifts = NULL, n_level = 1, unscaled_prob_cats, 
                                cats_num){
  # noise, res are p x n matrices!
  n_cats <- 0 
  for (j in 2:ncol(m)){
    mask <- which(m[,j] == 0) # to exclude not effective variables
    res_middle <- res
    res_middle[mask,] <- 0
    print(j)
    
    # check if needs one-hot encoding when cats are the parents
    if (!(j %in% idxy)){
      
      onehot <- intersect(idxy, which(m[,j] == 1))
      
      print(j+n_cats*(n_level-1))
      print(seeds[j+n_cats*(n_level-1)])
      res[j,] <- as.vector(mlp(t(res_middle), size, seed = seeds[j+n_cats*(n_level-1)], 
                               shift = shifts[j], onehot = onehot,
                               cat_level = 1)) + noise[j,] #  cat_level = 1 for numerical variables
      
    } else {
      # n_level - 1 more seeds for each coordinate

      print(j+n_cats*(n_level-1))
      print(j+(n_cats+1)*(n_level-1))
      print(seeds[(j+n_cats*(n_level-1)):(j+(n_cats+1)*(n_level-1))])
      res.pre <- as.matrix(mlp(t(res_middle), size, 
                               seed = seeds[(j+n_cats*(n_level-1)):(j+(n_cats+1)*(n_level-1))], 
                               shift = shifts[j], onehot = NULL, 
                               cat_level = n_level)) 
      
      res.pre <- res.pre + unscaled_prob_cats[which(idxy == j), , ]
      # COMPUTE ARGMIN
      res[j,] <- apply(res.pre,2,which.max)
      print(res[j,])
      n_cats <-  n_cats + 1 
    }
  }
  return(res)
}


### Simulating nonlinear data ###
sim.nonlinear <- function(n,p,ify = T, m = NULL, idxy = NULL, y_no_parent = F, 
                          simulator = 'mlp', seeds = NULL, size = 200, shifts = NULL,
                          n_level = 1, g.type = "ER", noise.type = "uniform", 
                          cats_num = 1){
  
  # n number of samples, p number of predictors
  
  # simulate the graph
  if (is.null(m)){
    if (g.type == "ER"){
      rdb <- rbinom(p*(p-1)/2, 1, 2/((p-1)))
      m <- diag(x=0, p)
      upperm = upper.tri(m, diag = F)
      m[upperm] = rdb
    } else if (g.type == "SF"){
      m.p <- 2
      g <- sample_pa(n = p, m = m.p, directed = T)
      m <- as.matrix(as_adjacency_matrix(g))
      m <- m[nrow(m):1,ncol(m):1]
    }
  } 
  
  # generate noises
  if (noise.type == "uniform"){
    noise <- matrix(runif(n*p, -3, 3), p, n) # simple noise 
  } else if (noise.type == "logis"){
    noise <- matrix(rlogis(n*p, 0, 2), p, n) # simple noise 
  }
  
  # sample position for Y
  if (ify){
  
    idxy <- sample(1:p, cats_num)
    # adjust the noise for Y, let's say NORMAL distribution
    # not used if n_level > 1
    for (i in idxy){
      noise[i,] <- rnorm(n, 0, 1)
    }
  } else if (is.null(idxy)){
    # THIS CORRESPONDS TO THE CASE WHEN THERE'S NO Y
    idxy <- -1
  } else {
    # for reproduction, idxy is specified
    idxy <- idxy
    noise[idxy,] <- rnorm(n, 0, 1)
  }
  print(idxy)
  res <- matrix(0, p, n)
  
  if (cats_num > 1){
    for (i in idxy){
      for (j in idxy[2:cats_num]){
        m[i, j] <- 0
        m[j, i] <- 0
      }
    }
  }
  
  child_cats <- list()
  parent_cats <- list()
  for (i in 1:cats_num){
    child_cats[[i]] <-which(m[idxy[i],] == 1)
    parent_cats[[i]] <- which(m[,idxy[i]] == 1)
  }
  
  if (n_level>1) {
    unscaled_prob_cats <- array(0, dim = c(cats_num, n_level, n))
    for (i in 1:cats_num){
      shift_cat <- sample(seq(-0.5,0.5,0.05), n_level, replace = T)
      print("shift_cat")
      print(shift_cat)
      # sample logis noise
      noise_cat <- matrix(rlogis(n*n_level), nrow = n_level)
      unscaled_prob_cats[i,,] <- sweep(noise_cat, 1, shift_cat, FUN = "+")
    }
    
    
  } else {
    unscaled_prob_cats <- NULL
    intercept.cats <- NULL
  }
  
  # adjust res and causal graph
  for (i in 1:cats_num){
    if (idxy[i] == 1){
      if (n_level == 1){
        res[1,] <- sign(noise[1,])
      } else {
        # sample shifts
        res[1,] <- apply(t(unscaled_prob_cats[i,,]),1,which.max) 
      }
      break
    }
    else{
      res[1,] <- noise[1,]
    }
  }
  
  
  # simulate the coefficients
  if (simulator == 'mlp'){
    if (is.null(seeds)){
      seeds <- sample(1:1000, p+(n_level-1)*cats_num)
    }
    print(paste("The current seeds are:"))
    print(seeds)
    res <- compute_samples_mlp(noise, idxy, res, seeds, m, size, shifts, 
                               n_level, unscaled_prob_cats, cats_num)
  } 
  
  samples <- t(res)
  
  print(m)
  
  ordered_idxy <- idxy[order(idxy)]
  for (i in 1:cats_num){
    current_id <- ordered_idxy[i] - (i - 1)
    # permutation to make Y the last variable
    m <- rbind(m[-idxy,], m[idxy,])
    m <- cbind(m[,-idxy], m[,idxy])
    
    noise <- rbind(noise[-current_id, ], noise[current_id, ])
    samples <- cbind(samples[, -current_id], samples[, current_id])
  }
  
  
  # generate Y
  return(list(graph = m, noise = noise, samples = samples, idxy = idxy))
  
}
