rm(list=ls()); gc()
library(rstan)
library(coda)
library(matrixStats)
library(data.table)
library(lme4)
library(MCMCglmm)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")

### We need to compute a mean phenotypic value for each line/population

final_export$temperature = (final_export$temperature-mean(final_export$temperature))/sd(final_export$temperature)
final_export$rel_humidity = (final_export$rel_humidity-mean(final_export$rel_humidity))/sd(final_export$rel_humidity)
final_export$logD = (final_export$logD-mean(final_export$logD))/sd(final_export$logD)

Lines_means=NULL
for(i in 1:nb_trait){
  temp_mod1 <- lmer(final_export[, paste0(vect_P_traits[i],"_pred_02")] ~ (temperature+rel_humidity+logD)^3 + pop_label - 1 + (1|date_str),data= final_export)
  Lines_means <- cbind(Lines_means,as.numeric(summary(temp_mod1)$coef[,1]))
}

row.names(Lines_means)=tstrsplit(names(summary(temp_mod1)$coef[,1]),"pop_label")[[2]]
Lines_means=Lines_means[4:(nrow(Lines_means)-4),]

mean_phen_values <- data.frame(line=rownames(Lines_means),
                               T12= Lines_means[,1],T13= Lines_means[,2],T21= Lines_means[,3],
                               T23= Lines_means[,4],T31= Lines_means[,5],T32= Lines_means[,6])

# Load the fertility table from Noble 2017
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
is_pop=!is.na(tstrsplit(mean_phen_values$line,'_')[[2]])
sav_Means_P_values <- as.numeric(colMeans(mean_phen_values[! is_pop,2:7]))

k=0
for( i in 9:14){
  k=k+1
  final_fertility[,i] <- final_fertility[,i] - sav_Means_P_values[k]
}

Y=final_fertility$fertility
Y_permuted <- NULL
for(i in 1:500) Y_permuted <- cbind(Y_permuted,sample(Y))
post_mod_gamma_Rdm=NULL
for(i in 1:500){
  print(i)
  Y_temp= Y_permuted[,i]
  model_MCMC_Rdm <- MCMCglmm(Y_temp ~ 1 + (T12+T13+T21+T23+T31+T32)^2 - (T12+T13+T21+T23+T31+T32) + I(T12^2) + I(T13^2) + I(T21^2) + I(T23^2) + I(T31^2) + I(T32^2), 
                             data = final_fertility, verbose = FALSE, nitt=150000, burnin=100000, thin=10)
  if(length(model_MCMC_Rdm)==20) post_mod_gamma_Rdm = cbind(post_mod_gamma_Rdm, posterior.mode(model_MCMC_Rdm$Sol)[2:22])
}

##############################################################################
################## Model with the non-permuted data ##########################
model_MCMC <- MCMCglmm(Y ~ 1 + (T12+T13+T21+T23+T31+T32)^2 - (T12+T13+T21+T23+T31+T32) + I(T12^2) + I(T13^2) + I(T21^2) + I(T23^2) + I(T31^2) + I(T32^2), 
                       data = final_fertility ,verbose = TRUE,nitt=50500000, burnin=500000,thin=100)

##############################################################################
##############################################################################

post_mod_gamma = posterior.mode(model_MCMC$Sol)[2:22]
temp_vect = post_mod_gamma[c(7:21,1:6)]
gamma <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                  temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                  temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                  temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                  temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                  temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

eigen(gamma)
proj_M = eigen(gamma)$vectors

extract_gamma_EV <- function(temp_vect){
  temp_vect = temp_vect[c(7:21,1:6)]
  gamma_rdm <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                        temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                        temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                        temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                        temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                        temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)
  
  # Rotation
  rot_gamma_rdm <- I(t(proj_M) %*% gamma_rdm %*% proj_M)
  # Selection gradient along the trait
  return(diag(rot_gamma_rdm))
}

rdm_EV_along_y=t(apply(post_mod_gamma_Rdm, 2, extract_gamma_EV))

plot(1:6, eigen(gamma)$values, 
     xlab=expression(paste("Canonical axis of selection (",y[i],")")),
     ylab=expression(paste("Strength of selection (eigenvalues ",lambda[i],")")),
     cex.axis=1.2,las=1,cex.lab=1.1,space=0,width=1,xlim=c(0.5,6.5),ylim=c(-13,13),bty="n",pch=8,type="n")
abline(h=0,lty=2)

all_rdm_EV_along_y_90 = t(apply(rdm_EV_along_y,2,function(x){ HPDinterval(as.mcmc(x),prob=.95) }))
all_rdm_EV_along_y_80 = t(apply(rdm_EV_along_y,2,function(x){ HPDinterval(as.mcmc(x),prob=.8) }))

arrows((1:6), all_rdm_EV_along_y_90[,1],(1:6), all_rdm_EV_along_y_90[,2],code=3,angle=90,length=.05,col="grey")
arrows((1:6), all_rdm_EV_along_y_80[,1],(1:6), all_rdm_EV_along_y_80[,2],code=3,angle=90,length=0,col="red",lwd=2)
points(c(1:6),colMedians(rdm_EV_along_y),pch=16,col="black")
points(1:6,eigen(gamma)$values,pch=8)
legend(0.3,14,c("Null dist. (median, 80%\n and 95% CI)","True est."),lwd=1,pch=c(16,8),bty="n",cex=.85)

####### NOW THE GAMMA ESTIMATES

### Intervals ###
post_dist_all <- HPDinterval(model_MCMC$Sol,prob=.95)[2:22,]
post_dist_all = post_dist_all[c(7:21,1:6),]

post_dist_all2 <- HPDinterval(model_MCMC$Sol,prob=.8)[2:22,]
post_dist_all2 = post_dist_all2[c(7:21,1:6),]

mode_post_Gamma <- posterior.mode(model_MCMC$Sol)[2:22]
mode_post_Gamma = mode_post_Gamma[c(7:21,1:6)]

##########################

plot(mode_post_Gamma,c(21:1),yaxt="n",bty="n",ylab="",xlim=c(-6,7),xlab="")
abline(v=0)
axis(side=2,at=21:1,labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                             "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                             "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                             "SF*SF","SB*SB","FS*FS","FB*FB","BS*BS","BF*BF"),las=1)

arrows(post_dist_all[,1],c(21:1),post_dist_all[,2],c(21:1),code=3,length=.05,angle=90)
arrows(post_dist_all2[,1],c(21:1),post_dist_all2[,2],c(21:1),code=3,length=0,angle=90,lwd=2,col="red")
points(mode_post_Gamma,c(21:1),pch=21,bg="black")

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")
