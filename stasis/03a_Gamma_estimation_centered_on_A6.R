rm(list=ls())
gc()
#load libraries
library(rstan)
library(coda)
library(matrixStats)
library(data.table)
library(lme4)

# Compute the Gamma matrix
# Includes the Figure 7 and S9


### THE STAN FUNCTION
getStanLinearModel = function() {
'
data {
  // number of observations
  int N;
  // response
  vector[N] y;
  // number of columns in the design matrix X
  int K;
  // design matrix X
  // should not include an intercept
  matrix [N, K] X;
  // keep responses
  int use_y_rep;
  int use_log_lik;
}
parameters {
  // regression coefficient vector
  real alpha;
  vector[K] beta;
  real sigma;
}
model {	
  // likelihood
  y ~ normal(alpha + X * beta, sigma);
}

'   
}

### END OF THE STAN FUNCTION

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")


### We need to compute a mean phenotypic value for each line/population

final_export$temperature= (final_export$temperature-mean(final_export$temperature))/sd(final_export$temperature)
final_export$rel_humidity= (final_export$rel_humidity-mean(final_export$rel_humidity))/sd(final_export$rel_humidity)
final_export$logD= (final_export$logD-mean(final_export$logD))/sd(final_export$logD)

Lines_means=NULL
for(i in 1:nb_trait){

temp_mod1 <- lmer(final_export[, paste0(vect_P_traits[i],"_pred_02")]~ (temperature+rel_humidity+logD)^3 +pop_label-1+(1|date_str),data= final_export)
Lines_means <- cbind(Lines_means,as.numeric(summary(temp_mod1)$coef[,1]))
}

row.names(Lines_means)=tstrsplit(names(summary(temp_mod1)$coef[,1]),"pop_label")[[2]]
Lines_means=Lines_means[4:(nrow(Lines_means)-4),]

mean_phen_values <- data.frame(line=rownames(Lines_means),
T12= Lines_means[,1],T13= Lines_means[,2],T21= Lines_means[,3],
T23= Lines_means[,4],T31= Lines_means[,5],T32= Lines_means[,6])

#Load the fertility table from Noble 2017
load("~/PATH/TO/DIR/fertility.rda")
	fertility <- subset(ecoefs,env=="NGM")

fertility_transition_traits <- subset(final_export,pop_label%in% fertility$line)
dim(fertility_transition_traits) # 556

fertility = subset(fertility,line%in% mean_phen_values$line)
fertility$line <- as.factor(as.character(fertility$line))
length(unique(fertility$line)) #232

final_fertility <- merge(fertility, mean_phen_values)
dim(final_fertility) # 232 OK

final_fertility$fertility <- ( exp(final_fertility$fertility)/mean(exp(final_fertility $fertility)))

## We will center using the mean A6140 value

sav_Means_P_values <- colMeans(mean_phen_values[,2:7])#colMeans(final_fertility[,9:14])
sav_Means_P_values =colMeans(sav_Means_P_values)
k=0
for( i in 9:14){
	k=k+1
	 final_fertility[,i] <- final_fertility[,i] - sav_Means_P_values[k]

}
attach(final_fertility)

X <- cbind(T12*T13,T12*T21,T12*T23,T12*T31,T12*T32,
T13*T21,T13*T23,T13*T31,T13*T32,T21*T23,T21*T31,T21*T32,T23*T31,T23*T32,T31*T32,
(final_fertility[,9:14])^2)

detach(final_fertility)

Y=final_fertility$fertility
hist(exp(Y))
mod1_data <- list(
  X = X,
  K = ncol(X),
  N = nrow(X),
  y = Y,
  use_y_rep = FALSE,
  use_log_lik = FALSE
)

stan_params = list(verbose = FALSE, chains = 1, warmup = 10000, iter = 15000, cores = 2, refresh = 100)

## Call the stan model

results_stan <- do.call(stan, c(
    list(model_code = getStanLinearModel(), 
         data = mod1_data,control=list(max_treedepth=15)), 
    stan_params))

check_treedepth(results_stan)
check_energy(results_stan)
print(results_stan)

post_beta <- As.mcmc.list(results_stan,pars=c("alpha","beta"))
post_dist_all <-HPDinterval(post_beta[[1]])

Y_axis_labels <- c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
"SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
"FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
"SF*SF","SB*SB","FS*FS","FB*FB","BS*BS","BF*BF")
names(results_stan)[2:22] <- Y_axis_labels

quartz()
#plot(results_stan,pars=c("alpha",Y_axis_labels))
plot(results_stan,pars=Y_axis_labels)

mode_estimates <- posterior.mode(post_beta[[1]])

temp_vect <- mode_estimates[2:22]
gamma <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

temp_vect = post_dist_all[2:22,1]

gamma_low <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

temp_vect = post_dist_all[2:22,2]

gamma_high <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)


write.table(matrix(paste0(round(1000* gamma)/1000
,"  [",round(1000* gamma_low)/1000,";",
round(1000*gamma_high)/1000,"]"),ncol=6),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file="~/PATH/TO/DIR/Cemee_Pop_WI/gamma.txt")

### We can produce a distribution of gammas from the stan model..
# Merge all 4 distributions

random_sampling <- matrix(floor(runif(21*5000,min=1,max=1001)),ncol=21)

eigen_rdm <- NULL
for(i in 1:5000){
temp_vect <- diag(post_beta[[1]][random_sampling[i,],2:22])

eigen_rdm=rbind(eigen_rdm,
eigen(matrix(c(temp_vect[16],temp_vect[1:5],
temp_vect[1],temp_vect[17],temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18],temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19],temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20],temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]),6,6))$values
)
}

# Figure 7
barplot(eigen_rdm_df$EV,xlab=expression(paste("Canonical axis of selection (",y[i],")")),ylab=expression(paste("Strength of selection (eigenvalues ",lambda[i],")")),cex.axis=1.2,las=1,cex.lab=1.1,space=0,width=1,xlim=c(0,7),ylim=c(-13,8.2))
vect_arrows =cbind(eigen_rdm_df$EV-eigen_rdm_df$sds,eigen_rdm_df$EV+eigen_rdm_df$sds)
arrows(seq(.5,5.5)+seq(0,0,by=.2),vect_arrows[,1],seq(.5,5.5)+seq(0,0,by=.2),vect_arrows[,2],code=3,angle=90,length=.05,lwd=1,col=grey(.5))
axis(side=1,at=seq(0.5,5.5,length=6),labels=1:6,pos=-14)
# Rotation matrix containing the eigen-vectors of gamma
proj_M <- eigen(gamma)$vectors

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

write.table(gamma,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file="~/PATH/TO/DIR/Cemee_Pop_WI/G_mat_for_simulations/gamma.txt")

