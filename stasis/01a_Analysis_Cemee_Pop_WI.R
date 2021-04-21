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

#Load populations data

final_merged =read.table("~/PATH/TO/DIR/Cemee_Pop_WI/Final_merged_data_NGM.txt",h=TRUE,sep="\t")

#We remove the last block which is very specific and will be analyzed separately
final_merged=subset(final_merged,!data_group_name=='B400')

vect_WI <- c("AB1","CB4507", "CB4852","CB4855" , "CB4856" ,"CB4858", "JU319","JU345" , "JU400" , "MY1","MY16","N2anc","PB306" ,"PX174","PX179" , "RC301")
vect_populations <- c("A6140", "CA150", "CA250", "CA350", "CA1100", "CA2100", "CA3100")

final_merged$population = as.factor((tstrsplit(final_merged$pop_label, split = "L", fixed = TRUE)[[1]]))
unique(final_merged$population)

#### We need to perform a correction for the user and the temperature/humidity
## Select the lines that have been performed in each location
# Blocks 301 to 309 are populations and not lines

table(unique(final_merged[!final_merged$data_group_name%in%paste0("B30",1:9),c("pop_label","location_label")])$location_label)
#Lisbon  Paris Paris2 
#   352    181    248

Lisbon=unique(subset(final_merged,location_label=="Lisbon")$pop_label)
Paris=unique(subset(final_merged,location_label=="Paris")$pop_label)
Paris2=unique(subset(final_merged,location_label=="Paris2")$pop_label)
all_loc=Lisbon[Lisbon%in%Paris & Lisbon%in%Paris2];length(all_loc)
#n=26

library(lme4)
all_loc=subset(final_merged,pop_label%in%all_loc)
all_loc$date_lab=as.character(substring(all_loc$date_str,1,4))

final_merged_predicted=final_merged
Exp1= (final_merged_predicted$location_label=="Paris")
Exp2= (final_merged_predicted $location_label=="Lisbon")
Exp3= (final_merged_predicted $location_label=="Paris2")

vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")

for(i in 1:length(vect_P_traits)){

vect_phen= all_loc[, vect_P_traits[i]]
all_loc_model=lmer(vect_phen ~pop_label+location_label+(1|date_str),data= all_loc)
coefs_fitted=as.numeric(coef(all_loc_model)[[1]][1,27:28])

final_merged_predicted = cbind( final_merged_predicted , final_merged_predicted[, vect_P_traits[i]] - 
									Exp1 * coefs_fitted[1] -
									Exp3 * coefs_fitted[2])
names(final_merged_predicted)[ncol(final_merged_predicted)]=paste0(vect_P_traits[i],"_pred_01")
}

final_merged_predicted$population = as.factor((tstrsplit(final_merged_predicted$pop_label, split = "L", fixed = TRUE)[[1]]))

## Second correction

Lisbon_2012=unique(subset(final_merged_predicted, substring(final_merged_predicted$date_str,1,4)=="2012")$pop_label)

Paris_all=unique(subset(final_merged_predicted, substring(final_merged_predicted$date_str,1,4)!="2012")$pop_label)

shared_2ndcor= Lisbon_2012[Lisbon_2012%in% Paris_all] 
length(shared_2ndcor) # n=64
shared_2ndcor=subset(final_merged_predicted,pop_label%in% shared_2ndcor)						
table(shared_2ndcor$location_label , substring(shared_2ndcor $date_str,1,4))

is_2012 = substring(final_merged_predicted$date_str,1,4)!="2012"	

for(i in 1:length(vect_P_traits)){

vect_phen= shared_2ndcor[,paste0(vect_P_traits[i],"_pred_01")]

shared_2ndcor_model=lmer(vect_phen ~ pop_label+(substring(date_str,1,4)=="2012") + (1|date_str),data= shared_2ndcor)
coefs_fitted = as.numeric(coef(shared_2ndcor_model)[[1]][1,65])

final_merged_predicted = cbind( final_merged_predicted , final_merged_predicted[, paste0(vect_P_traits[i],"_pred_01")] + is_2012 * coefs_fitted[1])
names(final_merged_predicted)[ncol(final_merged_predicted)]=paste0(vect_P_traits[i],"_pred_02")
}

### Corrections are done

#### FINAL EXPORT - SAME AS FOR MA LINES

final_export = subset(final_merged_predicted, select = c("data_id", "pop_label", "data_group_name", 
	"date_str", "time_str", "exper_duration_s", "logD", "temperature", "rel_humidity","env_label", vect_P_traits ,paste0(vect_P_traits,"_pred_01"),paste0(vect_P_traits,"_pred_02")))

final_export$population = as.factor((tstrsplit(final_export$pop_label, split = "L", fixed = TRUE)[[1]]))
final_export = final_export[, c(1, 2,29, 3:28)]


###### Another populations column
final_export$population2 <- as.character(final_export $population)
final_export$population2[final_export$population2%in%c("CA150","CA250","CA350")] <-  "CA50"
final_export$population2[final_export$population2%in%c("CA1100","CA2100","CA3100")] <-  "CA100"
final_export$population2[final_export$population2%in%c("GA150","GA250","GA450")] <-  "GA50"
final_export$population2[final_export$population2%in%vect_WI] <-  "WI"
final_export$population2[!final_export$population2=='WI' & final_export$data_group_name%in%c(paste0("B30",1:9),"B310")]<-  "pops"
final_export$population2 <- as.factor(final_export$population2)
table(final_export$population2)

write.table(final_export, file = "~/PATH/TO/DIR/Cemee_Pop_WI/Final_merged_Corrected_NGM.txt", 
	row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


rm(all_loc, all_loc_model, coefs_fitted,Exp1,Exp2,Exp3, final_merged,final_merged_predicted,i,is_2012,Lisbon,Lisbon_2012,Paris,Paris_all,Paris2,shared_2ndcor, shared_2ndcor_model,vect_phen)

data_for_G_estimation <- subset(final_export,population2%in%c("A6140","CA50","CA100"))

for(j in c('temperature',"rel_humidity","logD")){
data_for_G_estimation[,j] <- (data_for_G_estimation[,j]-mean(data_for_G_estimation[,j]))/sd(data_for_G_estimation[,j])
}



for(i in c(24:29)){
# The correction are performed for each population  individually
for(vect_pop_id in c("A6140","CA50","CA100")){

vec_phen <- data_for_G_estimation[data_for_G_estimation$population2== vect_pop_id,i]
sub_data <- data_for_G_estimation[data_for_G_estimation$population2== vect_pop_id,]
sub_data$pop_label=as.factor(as.character(sub_data$pop_label))
sub_data$date_str =as.factor(as.character(sub_data$date_str))
	model_Pop <- lmer(vec_phen~temperature+rel_humidity+logD+pop_label+(1|date_str),data= sub_data )
	env_effects <- summary(model_Pop)$coef[c(2,3,4)]


vect_temp <- data_for_G_estimation$temperature[data_for_G_estimation$population2== vect_pop_id]
vect_humid <- data_for_G_estimation$rel_humidity[data_for_G_estimation$population2== vect_pop_id]
vect_dens <- data_for_G_estimation$logD[data_for_G_estimation$population2== vect_pop_id]

data_for_G_estimation[data_for_G_estimation$population2== vect_pop_id,i]<- data_for_G_estimation[data_for_G_estimation$population2== vect_pop_id,i] - 
						env_effects[1] * vect_temp -
						env_effects[2] * vect_humid- 
						env_effects[3] * vect_dens		
}}
rm(env_effects,i,j,model_Pop,sub_data,vec_phen,vect_dens,vect_humid,vect_temp,vect_pop_id)

plot.acfs <- function(x) {
	n <- dim(x)[2]
	par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
	for (i in 1:n) {
		acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(-0.15, 0.15))
		grid()
	}
}

VCV_mat = NULL
nb_trait = length(vect_P_traits)
k=0

for (i in vect_populations) {
print(i)
	temp_final = subset(data_for_G_estimation,   population == i)

	phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12_pred_02, T13_pred_02, T21_pred_02, T23_pred_02, T31_pred_02, T32_pred_02)) ~  trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000)

	pdf(file = paste0("~/PATH/TO/DIR/Cemee_Pop_WI/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_", 
		i, ".pdf"), height = 20)
	par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
	plot(model_MCMC$Sol, auto.layout = F)
	dev.off()

	pdf(file = paste0("~/PATH/TO/DIR/Cemee_Pop_WI/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_autocorr_", 
		i, ".pdf"), height = 10)
	plot.acfs(model_MCMC$Sol)
	dev.off()

	post_dist = posterior.mode(model_MCMC$VCV)
	k=k+1
	VCV_mat[[k]]=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
}


###Create tables
k=0
for(vect_pop_id in c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")){
	k=k+1
write.table(matrix(paste0(round(1000*VCV_mat[[k]]$G1_mat)/1000
,"  [",round(1000*matrix(HPDinterval(VCV_mat[[k]]$VCV_Mat)[1:36,1],ncol=6))/1000,";",
round(1000*matrix(HPDinterval(VCV_mat[[k]]$VCV_Mat)[1:36,2],ncol=6))/1000,"]"),ncol=6),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("~/PATH/TO/DIR/Cemee_Pop_WI/G_mat_tables/", vect_pop_id,".txt"))

}

save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")
#

rm(list=ls())
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")
##### Statistics for the paper : first remove duplicates due to herm/male
is_male=(tstrsplit(final_export$data_id,"_",fixed=TRUE)[[4]]=="male")
is_male[is.na(is_male)]=FALSE
temp_df = subset(final_export,! is_male)
nrow(temp_df) # 1703 plates
length(unique(temp_df$date_str)) # 198 blocks
median(as.numeric(table(temp_df$date_str))) # median=8

length(unique(subset(final_export,population2=="A6140")$pop_label)) # 189
length(unique(subset(final_export,population2=="CA50")$pop_label)) # 155
length(unique(subset(final_export,population2=="CA50" & population=="CA150")$pop_label)) # 52
length(unique(subset(final_export,population2=="CA50" & population=="CA250")$pop_label)) # 52
length(unique(subset(final_export,population2=="CA50" & population=="CA350")$pop_label)) # 51

length(unique(subset(final_export,population2=="CA100")$pop_label)) # 172
length(unique(subset(final_export,population2=="CA100" & population=="CA1100")$pop_label)) # 51
length(unique(subset(final_export,population2=="CA100" & population=="CA2100")$pop_label)) # 53
length(unique(subset(final_export,population2=="CA100" & population=="CA3100")$pop_label)) # 68

temp_df=subset(final_export,population2%in%c("A6140","CA50","CA100"))
temp_df$pop_label=as.factor(as.character(temp_df$pop_label))
summary(as.numeric(table(temp_df$pop_label)))

length(unique(subset(final_export,population2=="GA50")$pop_label)) # 164
temp_df=subset(final_export,!population2%in%c("pops","WI"));temp_df$pop_label=as.factor(as.character(temp_df$pop_label))
summary(table(as.numeric(table(temp_df$pop_label))))
#  1   2   3   4   5   6 
# 30 540  55  28  26   1


### Populations and Wild isolates
length(unique(subset(final_export,population2!="GA50")$date_str))
length(unique(subset(final_export,population2=="WI")$date_str))

table(table(subset(final_export,population2=="WI")$pop_label))
table(table(subset(final_export,population2=="pops")$pop_label))

dim(unique(final_export[,c("population2","pop_label","date_str","time_str")]))
dim(final_export)

df_for_table_export=unique(final_export[,c("population2","pop_label","date_str","time_str")])
df_for_table_export$population2=as.character(df_for_table_export$population2)
df_for_table_export$population2[!df_for_table_export$population2%in%c("pops","WI")]="Inbred_line"
df_for_table_export$population2[df_for_table_export$population2=="pops"]="Population"
names(df_for_table_export)=c("ID1","ID2","Date","Time")
write.table(df_for_table_export,file="~/PATH/TO/DIR/Cemee_Pop_WI/Table_Assays.txt",col.names=TRUE,quote=FALSE,row.names=FALSE,sep="\t")