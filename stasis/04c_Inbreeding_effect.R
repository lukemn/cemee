rm(list=ls());gc()
library(MCMCglmm)
library(psych)
library(ggplot2)
library(dplyr)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(corrplot)
library(Rmisc)
library(nlme)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_with_gamma.RData")

### Another column
names(mean_phen_values)[1]<-'pop_label'
mean_phen_values=merge(mean_phen_values,unique(final_export[,c("pop_label","population",'population2')]))

mean_phen_values$pop_id <- mean_phen_values$pop_label
vec_temp_01 <- tstrsplit(mean_phen_values$pop_label,"_",fixed=TRUE)
vec_temp_02 <- tstrsplit(vec_temp_01[[1]],"L",fixed=TRUE)

vect_evol_pop_final <- tstrsplit(vec_temp_02[[1]],"noM",fixed=TRUE)[[1]]
vect_lines_or_pop_final <- vec_temp_01[[2]]

vect_lines_or_pop_final[is.na(vect_lines_or_pop_final) & as.numeric(vec_temp_02[[2]])==0] <- "WI"
vect_lines_or_pop_final[is.na(vect_lines_or_pop_final) & !as.numeric(vec_temp_02[[2]])==0] <- "cemee-lines"
mean_phen_values$pop_id = paste(vect_evol_pop_final, vect_lines_or_pop_final,sep="_")
## Remove "noMs"
remove_noMs <- tstrsplit(mean_phen_values$pop_label,"noM",fixed=TRUE)
mean_phen_values=mean_phen_values[is.na(remove_noMs[[2]]),]

table(mean_phen_values$population2)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# We need the confidence interval of the mean estimates for the populations only

vect_Cemee_pops=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")
Cemee_pops=subset(mean_phen_values,pop_id%in%paste0(vect_Cemee_pops,"_herm"))

names_model_coefs=tstrsplit(names(summary(temp_mod1)$coef[,1]),"pop_label")[[2]]
vect_names_lmer=names(summary(final_trait_model_list[[1]])$coef[,1])

mean_est_CI=list()
for(i in 1:7){
print(i)
mean_est_CI[[i]]=confint(final_trait_model_list[[i]],parm=vect_names_lmer[names_model_coefs%in%Cemee_pops$pop_label],nsim=100)
}

vect_traits=c("SF","SB","FS","FB","BS","BF")
mean_est_CI_df= NULL
for(i in 1:6){
mean_est_CI_df=rbind(mean_est_CI_df,data.frame(trait=vect_traits[i],pop_id=Cemee_pops$pop_label,mean_est_CI[[i]]))
}
names(mean_est_CI_df)=c("trait","pop_id","lwr","upr")


mean_est_CI_df$mean=NA
for(i in vect_Cemee_pops){
mean_est_CI_df[mean_est_CI_df$pop_id==paste0(i,"L0_herm"),]$mean=as.numeric(subset(mean_phen_values,pop_id==paste0(i,"_herm"))[,2:7])
	}

lines_predictions=NULL
for(i in 1:7){
	for(j in 1:6){
Cemee_lines=subset(mean_phen_values,pop_id==paste0(vect_Cemee_pops[i],"_cemee-lines"))
temp_df=data.frame(pop_id=paste0(vect_Cemee_pops[i],"_cemee-lines"),trait=names(Cemee_lines)[c(1+j)],unique(predict(lm(Cemee_lines[,(1+j)]~1),interval="confidence")))[1,]
lines_predictions=rbind(lines_predictions, temp_df)

}}
lines_predictions$pop_id=tstrsplit(lines_predictions$pop_id,"_")[[1]]
mean_est_CI_df$pop_id=tstrsplit(mean_est_CI_df $pop_id,"L0")[[1]]
names(mean_est_CI_df)[3:5]=paste0("pop_",names(mean_est_CI_df)[3:5])
names(lines_predictions)[3:5]=paste0("lines_",names(lines_predictions)[3:5])
lines_predictions$trait=levels(mean_est_CI_df $trait)[as.numeric((lines_predictions$trait))]

final_estimates=merge(mean_est_CI_df, lines_predictions)
head(final_estimates)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

v_col = c(grey(.5),"cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")
v_pch=c(21,22,24,3,4,23)

quartz(w=5,h=5)
par(mar=c(5,4,2,2))
plot(final_estimates $lines_fit~ final_estimates $pop_mean,las=1,bty="n",,ylim=c(-5,0),asp=1,xlab=("Mean population estimate"),ylab="Mean of all all lines estimates",type="n")
k=0
for(iter_pop in vect_Cemee_pops){
	temp_df=subset(final_estimates,pop_id==iter_pop)
	k=k+1
	arrows(temp_df$pop_mean,temp_df$lines_lwr,temp_df$pop_mean,temp_df$lines_upr,code=3,angle=90,length=.1,col=grey(.3))
arrows(temp_df$pop_lwr,temp_df$lines_fit,temp_df$pop_upr,temp_df$lines_fit,code=3,angle=90,length=.1,col=grey(.3))
}
abline(a=0,b=1)
k=0
for(iter_pop in vect_Cemee_pops){
	temp_df=subset(final_estimates,pop_id==iter_pop)

	k=k+1
	points(temp_df $lines_fit~temp_df$pop_mean,las=1,bty="n",col=v_col[k],bg=v_col[k],ylim=c(-5,0),pch=v_pch)
}
legend(-5,.3, c(vect_Cemee_pops,vect_traits,""),pch=c(rep(16,7), v_pch,0),col=c(v_col,rep("black",6),"white"),bty="n",ncol=2)

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_Inbreeding.RData")

