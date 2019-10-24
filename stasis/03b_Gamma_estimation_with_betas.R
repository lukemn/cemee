rm(list=ls())
gc()
#load libraries
library(rstan)
library(coda)
library(matrixStats)
library(data.table)
library(lme4)

# Produce figure S10

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

attach(final_fertility)
# The same as before with the Beta coefficients
X <- cbind(final_fertility[,9:14],T12*T13,T12*T21,T12*T23,T12*T31,T12*T32,
T13*T21,T13*T23,T13*T31,T13*T32,T21*T23,T21*T31,T21*T32,T23*T31,T23*T32,T31*T32,
(final_fertility[,9:14])^2)

detach(final_fertility)

Y=final_fertility$fertility

mod1_data <- list(
  X = X,
  K = ncol(X),
  N = nrow(X),
  y = Y,
  use_y_rep = FALSE,
  use_log_lik = FALSE
)

stan_params = list(verbose = FALSE, chains = 1, warmup = 3000, iter = 5000, cores = 2, refresh = 100)

## The stan model is here a function getStanLinearModel that you can find at the end of this file

results_stan <- do.call(stan, c(
    list(model_code = getStanLinearModel(), 
         data = mod1_data,control=list(max_treedepth=15)), 
    stan_params))

check_treedepth(results_stan)
check_energy(results_stan)
print(results_stan)

post_beta <- As.mcmc.list(results_stan,pars=c("alpha","beta"))
post_dist_all <-HPDinterval(post_beta[[1]])

Y_axis_labels <- c("SF","SB","FS","FB","BS","BF","SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
"SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
"FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
"SF*SF","SB*SB","FS*FS","FB*FB","BS*BS","BF*BF")
names(results_stan)[2:28] <- Y_axis_labels

quartz()

plot(results_stan,pars=Y_axis_labels)
median_estimates <- NULL
for(i in 1:28) median_estimates <- c(median_estimates,median(post_beta[[1]][,i]))

Y_pred_stan_ip <-  as.matrix(X)%*% median_estimates[2:28] + median_estimates[1]

temp_vect <- median_estimates[8:28]
gamma <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

temp_vect = post_dist_all[8:28,1]

gamma_low <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

temp_vect = post_dist_all[8:28,2]

gamma_high <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)


Z0 = - solve(gamma) %*% (post_dist_all[2:7])

post_dist_all[1] + 1/2* t(post_dist_all[2:7])%*%Z0
random_sampling <- matrix(floor(runif(21*5000,min=1,max=1001)),ncol=21)

eigen_rdm <- NULL
for(i in 1:5000){
temp_vect <-(post_beta[[1]][i,2:22])

eigen_rdm=rbind(eigen_rdm,
eigen(matrix(c(temp_vect[16],temp_vect[1:5],
temp_vect[1],temp_vect[17],temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18],temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19],temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20],temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]),6,6))$values
)
}

eigen_rdm_df <- data.frame(ndim=as.factor(1:6),EV=colMeans(eigen_rdm),EVmin=apply(eigen_rdm,2,function(x) sort(x)[round(.05*length(x))]),EVmax=apply(eigen_rdm,2,function(x) sort(x)[round(.95*length(x))]))


proj_M <- eigen(gamma)$vectors

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma_with_betas.RData")

