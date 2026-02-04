### generate synthetic datasets ###

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("data_simulationcode.R")
p_sets <- c(3, 5, 7, 10)
g.types <- c("ER", "SF")
noise.types <- c("uniform", "logis")
function.types <- c("linear", "MLP")
function.types <- c("linear")
dir.create(file.path("data/"), showWarnings = FALSE)
n_level <- 5 # choose the number of levels in categorical variables, e.g., 2, 3, 5


for (p in p_sets){
  for (g.type in g.types){
    for (noise.type in noise.types){
      for (function.type in function.types){
        for (rep in 1:20){
          set.seed(rep)
          if (function.type == "linear"){
            data.t <- sim.lingam2(1000, p, g.type = g.type, noise.type = noise.type,
                                  n_level = n_level, cats_num = floor(p/2))
          } else if (function.type == "MLP") {
            data.t <- sim.nonlinear(1000, p, g.type = g.type, noise.type = noise.type,
                                    n_level = n_level)
            # make sure the categorical variable has "n_level" levels in the outcome
            cats <- data.t$samples[,p]
            while (length(unique(cats)) != n_level){
              data.t <- sim.nonlinear(1000, p, g.type = g.type, noise.type = noise.type,
                                      n_level = n_level)
              cats <- data.t$samples[,p]
            }
          }
          # save data
          path <- paste0("data/sim_p",p,"_noise",noise.type,"_g",g.type,
                         "_function_",function.type, "_cats", floor(p/2), "_rep",rep,".rds")
          saveRDS(data.t,path)
        }
      }
    }
  }
}

