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
library(MuMIn)

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
  
  #is male frequency significant (unadjusted alpha=0.05)?
  if(anova(temp_mod1,temp_mod2,test="F")$Pr[2] < 0.05){
    retained_model[[i]] <- temp_mod1
    print(r.squaredGLMM(temp_mod1))
    lines(I(c(1:10)/10),predict(temp_mod1,data.frame(mf=I(c(1:10)/10)),re.form=NA),col="red")
    text(0.2,c(-2.5,-4,-2,-2,-1,1)[i],paste("p=",round(anova(temp_mod1,temp_mod2,test="F")$Pr[2],digits=3)))
    text(0.2,c(-2.5,-4,-2,-2,-1,1)[i],paste("p=",round(anova(temp_mod1,temp_mod2,test="F")$Pr[2],digits=3)))
    text(0.2,(c(-2.5,-4,-2,-2,-1,1)*1.25)[i],paste("r2=",round(r.squaredGLMM(temp_mod1)[1],digits=3)))
  }else{
    #no MF effect
    retained_model[[i]] <- temp_mod2
    print(r.squaredGLMM(temp_mod2))
    lines(I(c(1:10)/10),predict(temp_mod2,data.frame(is_pop=FALSE,mf=I(c(1:10)/10)),re.form=NA),col="black")
  }
}

vect_Cemee_pops=c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")

### END OF MALE EFFECTS

# Then we need to test the effect of inbreeding.
# We fit to the mean RIL estimate per population (not as random effect anymore)

data_to_test <- subset(final_inbreeding,population%in% vect_Cemee_pops)
data_to_test$pop_ispop <- paste(data_to_test$is_pop,data_to_test$population)

summary_results <- NULL
par(mfrow=c(2,3))
for(i in 1:6){
  final_model <- lmer(data_to_test[, paste0(vect_P_traits[i],"_pred_02")] ~ is_pop +(1|date_str)+(1|population/pop_label),data= data_to_test)
  final_model2 <- lmer(data_to_test[, paste0(vect_P_traits[i],"_pred_02")] ~1 +(1|date_str)+(1|population/pop_label),data= data_to_test)
  
  summary_results <- rbind(summary_results,
                           c(anova(final_model,final_model2)$Chisq[2],
                             anova(final_model,final_model2)$Pr[2]))
  qqnorm(resid(final_model))
  qqline(resid(final_model))
}

summary_results <- as.data.frame(summary_results)
names(summary_results) <- c("Chisq","P")
summary_results$traits <- vect_P_traits
summary_results <- summary_results[,c(3,1,2)]
summary(final_model)
