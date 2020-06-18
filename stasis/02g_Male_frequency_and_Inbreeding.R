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
library(lme4)

##  Effect of male frequency on the outbred population and Inbreeding effect
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")

male_freq <- read.table("~/PATH/TO/DIR/Cemee_Pop_WI/population_rep_male_freqs.txt",h=TRUE)
names(male_freq)[1] <- 'population'
head(male_freq)
# remove noM populations (no mf estimate)
male_freq  <- merge(male_freq[,c(1,2,5)] ,subset(data_populations,is.na(tstrsplit(data_populations$pop_label,"noM",fixed=TRUE)[[2]])))
dim(male_freq) # 102 37

## Subset the lines of interest

vect_keep <- (! final_export $population2%in%c("pops","WI") & substring(final_export$population,1,2)!="GA")

RILs_for_inbreeding <- subset(final_export,!population2%in%c("pops","WI") & substring(population,1,2)!="GA")
RILs_for_inbreeding$mf <- 0

final_inbreeding <- rbind(RILs_for_inbreeding, male_freq)

final_inbreeding$is_pop <- final_inbreeding$population2=="pops"

final_inbreeding_no_RILs <- subset(final_inbreeding,is_pop)
final_inbreeding_no_pops <- subset(final_inbreeding,!is_pop)

retained_model <- list()

par(mfrow=c(2,3))
for(i in 1:6){

vect_phenotypes= final_inbreeding_no_RILs[, paste0(vect_P_traits[i],"_pred_02")]
vect_phenotypes_all <- final_inbreeding[, paste0(vect_P_traits[i],"_pred_02")]
# Main model
temp_mod1 <- lmer(vect_phenotypes ~ (1|population)+(1|date_str)+mf,data= final_inbreeding_no_RILs)

# No male freq. effect
temp_mod2 <- lmer(vect_phenotypes ~ (1|population)+(1|date_str),data= final_inbreeding_no_RILs)


plot(vect_phenotypes_all ~ final_inbreeding$mf,pch=16,type="n",las=1,bty="n",ylab="Log transition rate (Hz)",xlab="Male frequency",main=c("SF","SB","FS","FB","BS","BF")[i])
points(vect_phenotypes ~ final_inbreeding_no_RILs$mf,pch=16,type="p")

#is mf significant ?
if(anova(temp_mod1,temp_mod2,test="F")$Pr[2] < 0.05){
retained_model[[i]] <- temp_mod1
lines(I(c(1:10)/10),predict(temp_mod1,data.frame(mf=I(c(1:10)/10)),re.form=NA),col="red")
text(0.2,c(-2.5,-4,-2,-2,-1,1)[i],paste("p=",round(anova(temp_mod1,temp_mod2,test="F")$Pr[2],digits=3)))
		
	
	}else{
	#no MF effect
	retained_model[[i]] <- temp_mod2
	lines(I(c(1:10)/10),predict(temp_mod2,data.frame(is_pop=FALSE,mf=I(c(1:10)/10)),re.form=NA),col="black")

	}
}


vect_Cemee_pops=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")

# Then we need to test the effect of inbreeding.
# We need to fit the mean RILs estimate per population (not as random effect anymore)


temp_df_traits = NULL

for(i in 1:nb_trait){

vect_phenotypes= final_inbreeding_no_RILs[, paste0(vect_P_traits[i],"_pred_02")]
temp_mod1 <- lmer(vect_phenotypes ~ population-1 +(1|date_str),data= final_inbreeding_no_RILs)

temp_df_traits <- rbind(temp_df_traits ,cbind(data.frame(trait=vect_P_traits[i],population= vect_Cemee_pops,mean_values =predict(temp_mod1,data.frame(population= vect_Cemee_pops),re.form=NA)),confint(temp_mod1,parm=paste0(rep("population",7), vect_Cemee_pops))[c(1,3,5,7,2,4,6),]))

}

vect_traits=c("SF","SB","FS","FB","BS","BF")
names(temp_df_traits)=c("trait","pop_id","pop_mean","pop_lwr","pop_upr")
temp_df_traits$trait <- vect_traits[as.numeric(temp_df_traits$trait)]

# Then we can perform the same with the RILs
temp_df_traits2 = NULL

for(i in 1:nb_trait){

vect_phenotypes= final_inbreeding_no_pops[, paste0(vect_P_traits[i],"_pred_02")]
temp_mod1 <- lmer(vect_phenotypes ~ population -1+(1|pop_label)+(1|date_str),data= final_inbreeding_no_pops)

temp_df_traits2 <- rbind(temp_df_traits2 ,cbind(data.frame(trait=vect_P_traits[i],population= vect_Cemee_pops,mean_values = predict(temp_mod1,data.frame(population= vect_Cemee_pops),re.form=NA)),confint(temp_mod1,parm=paste0(rep("population",7), vect_Cemee_pops))[c(1,3,5,7,2,4,6),]))

}

vect_traits=c("SF","SB","FS","FB","BS","BF")
names(temp_df_traits2)=c("trait","pop_id","lines_fit","lines_lwr","lines_upr")
temp_df_traits2$trait <- vect_traits[as.numeric(temp_df_traits2$trait)]

final_estimates <- merge(temp_df_traits, temp_df_traits2)


v_col = c(grey(.5),"cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")
v_pch=c(21,22,24,3,4,23)
v_pch2=c(21,22,24,3,4,23)[6:1]
quartz(w=5,h=5)

par(mar=c(5,4,2,2))
plot(final_estimates$lines_fit~ final_estimates $pop_mean,las=1,bty="n",ylim=c(-5,0),asp=1,xlab=("Mean population estimate"),ylab="Mean of all all lines estimates",type="n")
k=0
for(iter_pop in vect_Cemee_pops){
	temp_df=subset(final_estimates,pop_id==iter_pop)
	k=k+1
temp_col <- c(rgb(77, 77, 77, max = 255, alpha = 50),rgb(77, 77, 77, max = 255, alpha = 255))[as.numeric(sign(temp_df$lines_lwr-temp_df$pop_mean)==temp_df$lines_upr-temp_df$pop_mean)+1]


arrows(temp_df$pop_mean,temp_df$lines_lwr,temp_df$pop_mean,temp_df$lines_upr,code=3,angle=90,length=.1,col= temp_col)

temp_col <- c(rgb(77, 77, 77, max = 255, alpha = 50),rgb(77, 77, 77, max = 255, alpha = 255))[as.numeric(sign(temp_df$pop_lwr-temp_df$lines_fit)==sign(temp_df$pop_upr-temp_df$lines_fit))+1]

arrows(temp_df$pop_lwr,temp_df$lines_fit,temp_df$pop_upr,temp_df$lines_fit,code=3,angle=90,length=.1,col= temp_col)
}
abline(a=0,b=1)
k=0
for(iter_pop in vect_Cemee_pops){
	temp_df=subset(final_estimates,pop_id==iter_pop)
	temp_df <- temp_df[order(temp_df$trait),]
	k=k+1
	points(temp_df $lines_fit~temp_df$pop_mean,las=1,bty="n",col=v_col[k],bg=v_col[k],ylim=c(-5,0),pch=v_pch2)
}
legend(-5,.3, c(vect_Cemee_pops,vect_traits,""),pch=c(rep(16,7), v_pch,0),col=c(v_col,rep("black",6),"white"),bty="n",ncol=2)

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_Inbreeding.RData")
