rm(list=ls())
gc()
library(MCMCglmm)
library(stats)
library(parallel)
dist_list <- NULL

source('Simulations/Model_functions.R', chdir = TRUE)

final_param_df <- data.frame(nb_loci_scaling_factor=log(rep(seq(exp(.005), exp(.3), length.out=20),1)),lambda=3,nb_loci=100)

final_param_list=list()
for(i in 1:nrow(final_param_df)) final_param_list[[i]] <- final_param_df[i,]

no_cores <- 3
clust <- makeCluster(no_cores)

List_pop_mut_array <-parLapply(clust, final_param_list , Produce_Mut_Gen_arrays_for_G_simuls)
list_G_simulated <- parLapply(clust, List_pop_mut_array, produce_G_from_phen_array)
stopCluster(clust)

	A6140_mat <- read.table("~/PATH/TO/DIR/Cemee_Pop_WI/G_mat_for_simulations/A6140.txt", 
		sep = "\t")
	G_init_mat <- as.matrix(A6140_mat/2)
dist_list=NULL
for(i in 1:length(list_G_simulated)) dist_list <- rbind(dist_list,unlist(compute_G_distances(G_init_mat, list_G_simulated[[i]]$G1_mat/2)))
dist_df <- as.data.frame(dist_list)
names(dist_df) <- c("angle","asp","summed_Val")

dist_df$nb_loci_scaling_factor <- final_param_df$nb_loci_scaling_factor

quartz()
plot(dist_df$nb_loci_scaling_factor,dist_df$summed_Val,pch=16,cex=.6,col="black")

temp_results <-tapply(dist_df $summed_Val, dist_df$nb_loci_scaling_factor,mean)
points(as.numeric(names(temp_results)),as.numeric(temp_results),col="red",type="l",lwd=2)

temp_mod <- lm(summed_Val ~ nb_loci_scaling_factor,data= dist_df[6:9,])
abline(temp_mod,lwd=2,col="green")
abline(h=0);abline(v=-coef(temp_mod)[1]/coef(temp_mod)[2])



summed_diff=NULL
for(i in 1:length(list_G_simulated)){
summed_diff=c(summed_diff,sum(abs(G_init_mat - list_G_simulated[[i]]$G1_mat/2)))
}
plot(summed_diff~dist_df$nb_loci_scaling_factor)

first_estimate <- (-coef(temp_mod)[1]/coef(temp_mod)[2])
#0.1067742

save(list=ls(),file='~/PATH/TO/DIR/Simulations/Final_est_G_scaling_v2.RData')

### We have a first estimate, we want to refine

final_param_df <- data.frame(nb_loci_scaling_factor=rep(seq(.09, .13, length.out=100),1), lambda=3, nb_loci=100)

final_param_list=list()

for(i in 1:nrow(final_param_df)) final_param_list[[i]] <- final_param_df[i,]
no_cores <- 3

clust <- makeCluster(no_cores)
List_pop_mut_array <-parLapply(clust, final_param_list , Produce_Mut_Gen_arrays_for_G_simuls)
list_G_simulated2 <- parLapply(clust, List_pop_mut_array, produce_G_from_phen_array)
stopCluster(clust)

list_G_simulated= list_G_simulated2

dist_list=NULL
for(i in 1:length(list_G_simulated)) dist_list <- rbind(dist_list,unlist(compute_G_distances(G_init_mat, list_G_simulated[[i]]$G1_mat/2)))
dist_df <- as.data.frame(dist_list)
names(dist_df) <- c("angle","asp","summed_Val")

dist_df$nb_loci_scaling_factor <- final_param_df$nb_loci_scaling_factor

quartz()
plot(dist_df$nb_loci_scaling_factor,dist_df$summed_Val,pch=16,cex=.6,col="black")

temp_results <-tapply(dist_df $summed_Val, dist_df$nb_loci_scaling_factor,mean)
points(as.numeric(names(temp_results)),as.numeric(temp_results),col="red",type="l",lwd=2)

temp_mod <- lm(summed_Val ~ nb_loci_scaling_factor,data= subset(dist_df, nb_loci_scaling_factor>.09  & nb_loci_scaling_factor<.11 ))
abline(temp_mod,lwd=2,col="green")
abline(h=0);abline(v=-coef(temp_mod)[1]/coef(temp_mod)[2])
final_estimate <- (-coef(temp_mod)[1]/coef(temp_mod)[2])

summed_diff=NULL
summed_diff_diag=NULL
for(i in 1:length(list_G_simulated)){
summed_diff_diag =c(summed_diff_diag,sum(abs(diag(G_init_mat) - diag(list_G_simulated[[i]]$G1_mat)/2)))
summed_diff=c(summed_diff,sum(abs((G_init_mat) -(list_G_simulated[[i]]$G1_mat)/2)))
}
par(mfrow=c(1,2))
plot(summed_diff~dist_df$nb_loci_scaling_factor)
temp_mod2 <- lm(summed_diff ~ dist_df$nb_loci_scaling_factor+I(dist_df$nb_loci_scaling_factor^2))
points(predict(temp_mod2, dist_df) ~ dist_df$nb_loci_scaling_factor,col="green")
abline(h=0);abline(v=-coef(temp_mod)[1]/coef(temp_mod)[2])
abline(v=c(-summary(temp_mod2)$coef[,1][2]/summary(temp_mod2)$coef[,1][3]/2),col='green')

plot(summed_diff_diag ~dist_df$nb_loci_scaling_factor)
temp_mod3 <- lm(summed_diff_diag ~ dist_df$nb_loci_scaling_factor+I(dist_df$nb_loci_scaling_factor^2))
points(predict(temp_mod3, dist_df) ~ dist_df$nb_loci_scaling_factor,col="green")
abline(h=0);abline(v=-coef(temp_mod)[1]/coef(temp_mod)[2])
abline(v=c(-summary(temp_mod2)$coef[,1][2]/summary(temp_mod2)$coef[,1][3]/2),col='green')
abline(v=c(-summary(temp_mod3)$coef[,1][2]/summary(temp_mod3)$coef[,1][3]/2),col='red')

final_estimate=c(-summary(temp_mod2)$coef[,1][2]/summary(temp_mod2)$coef[,1][3]/2)
#  0.1086921
