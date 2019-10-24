rm(list = ls())
gc()

library(MCMCglmm)
library(psych)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(dplyr)


final_merged=read.table("~/PATH/TO/DIR/MA_lines/Final_merged_data_MA_Lines.txt",sep="\t",h=TRUE)
head(final_merged)

#Column IDs:
#exper_data_id and data_id are unique identifier of each movie
#pop_label is an unique identifier of each line
#data_group_name is an unique identifier of the block
#data_str is the date of the phenotyping in the format YYYYMMDD
#time_str is the time of the phenotyping in the format HHMMSS
#exper_duration_s is the length of each movie
#trackL is the mean track length
#n_tracks is the number of analyzed tracks
#mean_nb_tracks is the number of tracks normalized by the track length and the movie duration:
# trackL * n_tracks /exper_duration_s
#logD is the log10 of mean_nb_tracks
# temperature is the temperature at the beginning of each movie
# rel_humidity is the relative humidity at the beginning of each movie
# location_label is a character string that is not usefull for this data set (all values equal)
# env_label is a character string that is not usefull for this data set (all values equal)

# The final 6 columns are the 6 transition rates
vect_P_traits = c("T12", "T13", "T21", "T23", "T31","T32")
##  1,2,3 refers to the states Still, Forward and Backward so that "T12"
## stand for the transition rate from Still to Forward

nb_trait = length(vect_P_traits)

# Ancestral line of each MA line: either N2 or PB (short for PB306)
final_merged$anc_line=substring(final_merged$pop_label,1,2)

# A function to plot the autocorrelation from the MCMCglmm results
plot.acfs <- function(x) {
	n <- dim(x)[2]
	par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
	for (i in 1:n) {
		acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(0, 0.3))
		grid()
	}
}

VCV_mat = NULL

## We will now correct the block effects using the Ancestor data

df_ancestors <- subset(final_merged,pop_label%in%c("N2ancL0","PBancL0"))
block_id_for_corrections <- as.numeric(as.factor(final_merged $data_group_name))

for(i in 1: nb_trait){

#Extract the trait values for the ancestor lines
trait_values <- df_ancestors[,vect_P_traits[i]]
# Estimate the mean block effect
block_effect_model <- lm(trait_values ~ df_ancestors$data_group_name-1+df_ancestors$pop_label )

#correct the initial values
coefs_blocks=as.numeric(coef(block_effect_model)[1:16])[block_id_for_corrections]

#add the correct values to the initial data frame
final_merged <- cbind(final_merged, final_merged[,vect_P_traits[i]]+ coefs_blocks)
#name each corrected trait value so that T12 becomes T12_pred
names(final_merged)[ncol(final_merged)]=paste0(vect_P_traits[i],"_pred")

}

# Split the data set by ancestral lines and only retain the MA lines
N2lines= subset(final_merged, anc_line =="N2" & pop_label!="N2ancL0")
PBlines= subset(final_merged, anc_line =="PB" & pop_label!="PBancL0")

#Update the phenotypic trait vector
vect_P_traits=paste0(vect_P_traits,"_pred")

k=0
## We will now estimate the two M matrix using MCMCglmm for each set of MA lines separately
for (i in c("N2", "PB")) {

# First N2 and then PB306
	if(i=="N2") temp_df = N2lines
	if(i=="PB") temp_df = PBlines

# Transform the block column from a factor to a character string format
	temp_df$data_group_name = as.character(temp_df$data_group_name)

	
	phen.var = diag(nb_trait) * diag(var(subset(temp_df, select = vect_P_traits)))

	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12_pred, T13_pred, T21_pred, T23_pred, T31_pred, T32_pred )) ~ 
		(logD +rel_humidity+ temperature) + trait - 1, random = ~us(trait):pop_label + us(trait):data_group_name,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = temp_df, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000)

	pdf(file = paste0("~/PATH/TO/DIR/MA_lines/auto_corr_plots_MCMCglmm/M_for_Model_pdf_MCMC_", 
		i, ".pdf"), height = 20)
	par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
	plot(model_MCMC$Sol, auto.layout = F)
	dev.off()

	pdf(file = paste0("~/PATH/TO/DIR/MA_lines/auto_corr_plots_MCMCglmm/M_for_Model_pdf_MCMC_autocorr_", 
		i, ".pdf"), height = 10)
	plot.acfs(model_MCMC$Sol)
	dev.off()

	post_dist = posterior.mode(model_MCMC$VCV)
	k=k+1
	VCV_mat[[k]]=list(Population = i, N_measurement = nrow(temp_df), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)


}

dimnames(VCV_mat[[1]]$G1_mat) = list(vect_P_traits, vect_P_traits)
dimnames(VCV_mat[[2]]$G1_mat) = list(vect_P_traits, vect_P_traits)

write.table(VCV_mat[[1]]$G1_mat,file="~/PATH/TO/DIR/MA_lines/VM_estimates_N2.txt",sep="\t",quote=FALSE)
write.table(VCV_mat[[2]]$G1_mat,file="~/PATH/TO/DIR/MA_lines/VM_estimates_PB.txt",sep="\t",quote=FALSE)

save(list=ls(),file="~/PATH/TO/DIR/MA_lines/M_matrices_estimates.RData")

N2_mat <- VCV_mat[[1]]
PB_mat <- VCV_mat[[2]]

save(N2_mat ,file="~/PATH/TO/DIR/MA_lines/VCV_mat_N2.RData")
save(PB_mat ,file="~/PATH/TO/DIR/MA_lines/VCV_mat_PB.RData")

# Statistics for the manuscript

length(unique(N2lines$pop_label))
length(unique(PBlines$pop_label))

#### 

load("~/PATH/TO/DIR/MA_lines/M_matrices_estimates.RData")
vect_tnames=c("SF","SB","FS","FB","BS","BF")
df_for_barplot=data.frame(trait=final_merged[,23],anc_line= final_merged$anc_line,pop_label= final_merged$pop_label,tname= vect_tnames[1])
k=1
for(i in 24:28){
	k=k+1
	df_for_barplot=rbind(df_for_barplot,
data.frame(trait=final_merged[,i],anc_line= final_merged$anc_line,pop_label= final_merged$pop_label,tname= vect_tnames[k]))
}
par(mfrow=c(1,1),bty="n",las=1,xaxt="n")
boxplot(df_for_barplot$trait~ df_for_barplot$anc_line+ df_for_barplot$tname,col=rep(c("darkgreen","magenta")))
par(xaxt="s")
axis(side=1,at=c(1.5,3.5,5.5,7.5,9.5,11.5),labels=vect_tnames)






