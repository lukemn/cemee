rm(list=ls());gc()
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(data.table)
library(matrixStats)
library(boot)
library(Rmisc)
library(nlme)
library(parallel)
library(dae)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")


#Export phenotypic data for genomic analysis

run_parallel_MCMC <- function(list_param){
  i=list_param[[1]]
  temp_final=list_param[[2]]
  nb_trait=6
  vect_P_traits <- c("T12", "T13", "T21", "T23", "T31", "T32")
  
  temp_final$pop_label <- sample(temp_final$pop_label)
  temp_final$date_str <- sample(temp_final$date_str)	
  
  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), 
                             G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12_pred_02, T13_pred_02, T21_pred_02, T23_pred_02, T31_pred_02, T32_pred_02)) ~ trait-1, 
                         random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), 
                         data = temp_final, 
                         prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000)
  
  post_dist = posterior.mode(model_MCMC$VCV)
  
  VCV_mat_temp=list(Population = i, N_measurement = nrow(temp_final), 
                    G1_mat = matrix(post_dist[1:nb_trait^2], nb_trait, nb_trait), 
                    G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                    R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
  return(VCV_mat_temp)
  
}

clust <- makeCluster(25)

for(k in vect_populations){
  param_list=list()
  for(i in 1:1000) param_list[[i]] <- list(i=k, temp_final = subset(data_for_G_estimation, population == k))
  
  List_output <- parLapply(clust, param_list , run_parallel_MCMC)
  save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_.RData"))
  rm(List_output); gc()
}


## Subset of 50 A6140 RILs
uniq_A6140 <- as.character(unique(subset(data_for_G_estimation,population== vect_populations[1])$pop_label)) # 189

param_list=list()
for(i in 1:1000) param_list[[i]] <- list(i="A6140",temp_final = subset(data_for_G_estimation,   pop_label %in% sample(uniq_A6140,50)))

List_output <-parLapply(clust, param_list , run_parallel_MCMC)
rm(param_list)
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_subset.RData"))
rm(List_output);gc()
stopCluster(clust)
