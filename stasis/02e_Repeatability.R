rm(list=ls())
gc()
#load libraries
library(rstan)
library(coda)
library(matrixStats)
library(data.table)
library(lme4)
library(rptR)

# Generate the analysis and the plot of the repeatabilities (S7 and S10)

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")

### We need to compute a mean phenotypic value for each line/population

######### Inbred lines
rptr_mod_list=list()
k=0
for(j in levels(final_export$population2)[c(1,2,3)]){
temp_data=subset(data_for_G_estimation,population2==j)

	for(i in 1:nb_trait){
	k=k+1
	print(k)
	temp_data$vect_x=temp_data[, paste0(vect_P_traits[i],"_pred_02")]
	rptr_mod_list[[k]] <- rpt(vect_x~ (temperature+rel_humidity+logD) + (1|pop_label)+(1|date_str),datatype = "Gaussian",grname="pop_label",data= temp_data,parallel=TRUE,ncores=2)
}
}

######### Add the populations
load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData")
k=length(rptr_mod_list)
temp_data= subset(data_populations,population2=="pops" & !substring(pop_label,1,2)%in%c("CD","CM"))

	for(i in 1:nb_trait){
	k=k+1
	print(k)
	temp_data $vect_x=temp_data[, paste0(vect_P_traits[i],"_pred_02")]
	rptr_mod_list[[k]] <- rpt(vect_x ~ (temperature+rel_humidity+logD)+(1|population)+(1|data_group_name),datatype = "Gaussian",grname="population",data= temp_data,parallel=TRUE,ncores=2)
	}

### The WI
temp_data= subset(data_populations,population2=="WI")
	for(i in 1:nb_trait){
	k=k+1
	print(k)
	temp_data $vect_x=temp_data[, paste0(vect_P_traits[i],"_pred_02")]
	rptr_mod_list[[k]] <- rpt(vect_x ~ (temperature+rel_humidity+logD)+(1|population)+(1|data_group_name),datatype = "Gaussian",grname="population",data= temp_data,parallel=TRUE,ncores=2)
	}

## Create a DF
rptr_df=NULL
k=0
for(j in c(levels(final_export$population2)[c(1,2,3)],"pops","WI")){
	for(i in 1:nb_trait){
		k=k+1
		rptr_df=rbind(rptr_df,data.frame(pop=j,trait=vect_P_traits[i],R=as.numeric(rptr_mod_list[[k]]$R),rptr_mod_list[[k]]$CI_emp))
		}
}

names(rptr_df)=c("pop","trait","R","lwr","upr")
rptr_df$pop2=c(1,3,2,4,5)[as.numeric(as.factor(rptr_df$pop))]

#### And for the MA lines

load("~/PATH/TO/DIR/MA_lines/M_matrices_estimates.RData")

MA_rptr_df=NULL
for(j in 1:2){
	if(j==1) temp_data=N2lines
	if(j==2) temp_data=PBlines
	
	for(i in 1:nb_trait){
			print(i)
	temp_data$vect_x=temp_data[,vect_P_traits[i]]
		rptr_mod_temp =rpt(vect_x ~ (logD + rel_humidity +temperature) + (1|pop_label),datatype = "Gaussian",grname="pop_label",data= temp_data,parallel=TRUE,ncores=1)
		MA_rptr_df =rbind(MA_rptr_df,
		data.frame(pop=c("N2","PB306")[j],trait=vect_P_traits[i],rptr_mod_temp$R,rptr_mod_temp$CI_emp))
		}
}
names(MA_rptr_df)=c("pop","trait","R","lwr","upr")
MA_rptr_df$pop2=6
MA_rptr_df$pop2[MA_rptr_df$pop=="N2"]=7
MA_rptr_df$trait <- substring(MA_rptr_df$trait,1,3)

##########################################################
######### Figure B5 - RILS and MA lines ##################
##########################################################
rptr_df_plot1=rbind(subset(rptr_df,pop2%in%c(1:3)), MA_rptr_df)
par(mar=c(5,4,1,1))
vect_x= rep(1:5) + rep(c(0,8,16,24,32,40),each=5)
plot(rptr_df_plot1$R[order(rptr_df_plot1$trait,rptr_df_plot1$pop2)]~ vect_x ,ylab="Repeatability",bty="n",las=1,xaxt="n",ylim=c(0,1.1),xlab="",type="n")
arrows(vect_x, rptr_df_plot1$lwr[order(rptr_df_plot1$trait, rptr_df_plot1 $pop2)],vect_x, rptr_df_plot1$upr[order(rptr_df_plot1$trait, rptr_df_plot1 $pop2)],code=3,angle=90, length =.05)
points(rptr_df_plot1$R[order(rptr_df_plot1$trait, rptr_df_plot1$pop2)]~ vect_x,pch=21,bg= c(grey(.5),"cornflowerblue","firebrick3","white","darkgreen","magenta","darkgreen")[rptr_df_plot1$pop2[order(rptr_df_plot1$trait, rptr_df_plot1$pop2)]])
axis(side=1,at=c(3,11,19,27,35,43),labels=c("SF","SB","FS","FB","BS","BF"))

legend(2,1,c("Ancestor (A6140)","Generation 50 (CA[1-3]50)","Generation 100 (CA[1-3]100)","MA lines - PB306","MA lines - N2"),bty="n",pch=16, col=c(grey(.5),"cornflowerblue","firebrick3","magenta","darkgreen"),ncol=2)


#################################################################
######### Repeatabilities of the populations per window #########
#################################################################



pop_gen_df <- subset(data_populations,population2=="pops" & !substring(pop_label,1,2)%in%c("CD","CM"))
ref_gen <- data.frame(population=as.character(unique(pop_gen_df$population)),gen=c(0,rep(c(100,10,30,60),3),c(100,10,140,30,60),c(100,10,140,60),c(100,10,140,60),c(100,10,36,50,5,68),c(100,10,36,50,5,68),c(100,10,36,50,5,68),c(32,66,32,66,32,-5)))

ref_gen[27:49,]$gen <- ref_gen[27:49,]$gen + 140
pop_gen_df <- merge(pop_gen_df, ref_gen)

pop_gen_df$window <- 1
pop_gen_df$window[pop_gen_df$gen>=60] <- 2
pop_gen_df$window[pop_gen_df$gen==140] <- 3
pop_gen_df$window[pop_gen_df$gen>140] <- 4
pop_gen_df$window[pop_gen_df$gen>190] <- 5

k=0
rptr_mod_list_pop <- list()
for(j in unique(pop_gen_df$window)){
	temp_data <- subset(pop_gen_df,window==j)
	for(i in 1:nb_trait){
	k=k+1
	print(k)
	temp_data$vect_x=temp_data[, paste0(vect_P_traits[i],"_02")]
	rptr_mod_list_pop[[k]] <- rpt(vect_x ~ (temperature+rel_humidity+logD)+(1|population)+(1|data_group_name),datatype = "Gaussian",grname="population",data= temp_data,parallel=TRUE,ncores=2)
	}
}
# And the still frequency
k <- length(rptr_mod_list_pop)
for(j in unique(pop_gen_df$window)){
		temp_data <- subset(pop_gen_df,window==j)
		k=k+1
		temp_data$vect_x=temp_data[, "mean_still"]
	rptr_mod_list_pop[[k]] <- rpt(vect_x ~ (temperature+rel_humidity+logD)+(1|population)+(1|data_group_name),datatype = "Gaussian",grname="population",data= temp_data,parallel=TRUE,ncores=2)
}


rptr_df_pops=NULL
k=0
for(j in unique(pop_gen_df$window)){
	for(i in 1:nb_trait){
		k=k+1
		rptr_df_pops <-rbind(rptr_df_pops,data.frame(window=j,trait=vect_P_traits[i],R=as.numeric(rptr_mod_list_pop[[k]]$R), rptr_mod_list_pop[[k]]$CI_emp))
		}
}



names(rptr_df_pops)=c("window","trait","R","lwr","upr")
vectx <- c(30,90,140,190,240)
par(mfrow=c(2,3))
for(i in 1:6){
	temp <- subset(rptr_df_pops,trait==vect_P_traits[i])
	temp2 <- subset(rptr_df,trait==substring(vect_P_traits,1,3)[i] & pop=="WI")
plot(c(-33,vectx[temp$window]), c(temp2$R,temp $R),type="b",ylim=c(0,1),bty="n",las=1,pch=21,xlim=c(-50,260),bg=rep(c("darkgreen","grey"),c(1,5)),ylab="Repeatability",xlab="Generations")
arrows(vectx[temp$window], temp$lwr, vectx[temp$window], temp$upr,code=3,angle=90, length =.05)
arrows(-33, temp2$lwr, -33, temp2$upr,code=3,angle=90, length =.05)

}


temp_data= subset(data_populations,population2=="WI")
temp_data$vect_x=temp_data[, "mean_still"]
	rptr_WI_mean_still <- rpt(vect_x ~ (temperature+rel_humidity+logD)+(1|population)+(1|data_group_name),datatype = "Gaussian",grname="population",data= temp_data,parallel=TRUE,ncores=2)

rptr_df_Mean_Still=NULL

rptr_df_Mean_Still <-rbind(rptr_df_Mean_Still,data.frame(window=0,trait="mean_still",R=as.numeric(rptr_WI_mean_still$R), rptr_WI_mean_still$CI_emp))

# Figure 2D
k=30
for(j in unique(pop_gen_df$window)){
		k=k+1
		rptr_df_Mean_Still <-rbind(rptr_df_Mean_Still,data.frame(window=j,trait="mean_still",R=as.numeric(rptr_mod_list_pop[[k]]$R), rptr_mod_list_pop[[k]]$CI_emp))
		}
names(rptr_df_Mean_Still)=c("window","trait","R","lwr","upr")
par(mfrow=c(1,1))
rptr_df_Mean_Still$vect_x <- c(-33,30,90,140,190,240)
plot(rptr_df_Mean_Still$vect_x, rptr_df_Mean_Still$R,type="b",ylim=c(0,1),bty="n",las=1,pch=21,xlim=c(-50,260),bg=rep(c("darkgreen","grey"),c(1,5)),ylab="Repeatability",xlab="Generations")

arrows(rptr_df_Mean_Still$vect_x, rptr_df_Mean_Still$lwr, rptr_df_Mean_Still$vect_x, rptr_df_Mean_Still$upr,code=3,angle=90, length =.05)


save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Repeatability.RData")
