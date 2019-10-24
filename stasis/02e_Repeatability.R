rm(list=ls())
gc()
#load libraries
library(rstan)
library(coda)
library(matrixStats)
library(data.table)
library(lme4)
library("rptR")

# Generate the analysis and the plot of the repeatabilities

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
temp_data= subset(data_populations,population2=="pops")

	for(i in 1:nb_trait){
	k=k+1
	print(k)
	temp_data $vect_x=temp_data[, paste0(vect_P_traits[i],"_pred_02")]
	rptr_mod_list[[k]] <- rpt(vect_x ~ (temperature+rel_humidity+logD)+(1|population)+(1|data_group_name),datatype = "Gaussian",grname="population",data= temp_data,parallel=TRUE,ncores=2)
	}


## Create a DF
rptr_df=NULL
k=0
for(j in c(levels(final_export$population2)[c(1,2,3)],"pops")){
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

#### And now the plot

v_col = c(grey(.5),"cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

quartz(h=4,w=9)
layout(matrix(c(1,2),1,2),w=c(6,4))
par(mar=c(5,4,1,1))
vect_x= rep(1:4,6) + rep(c(0,8,16,24,32,40),each=4)
plot(rptr_df$R[order(rptr_df$trait)]~ vect_x ,ylab="Repeatabillity",bty="n",las=1,xaxt="n",ylim=c(0,1.1),xlab="",type="n")
arrows(vect_x,rptr_df$lwr[order(rptr_df$trait,rptr_df$pop2)],vect_x,rptr_df$upr[order(rptr_df$trait,rptr_df$pop2)],code=3,angle=90, length =.05)
points(rptr_df$R[order(rptr_df$trait,rptr_df$pop2)]~ vect_x,pch=21,bg= c(grey(.5),"cornflowerblue","firebrick3","white")[rptr_df$pop2[order(rptr_df$trait,rptr_df$pop2)]])
axis(side=1,at=c(4,12,20,28,36,44),labels=c("SF","SB","FS","FB","BS","BF"))

legend(2,1.2,c("Ancestor (A6140)","Generation 50 (CA[1-3]50)","Generation 100 (CA[1-3]100)"),bty="n",pch=c(16,16,16), col=c(grey(.5),"cornflowerblue","firebrick3"),ncol=1)
legend(27,1.2,"Populations",bty="n",pch=1,)


par(mar=c(5,0,1,1))
vect_x= rep(1:2,6) + rep(c(0,8,16,24,32,40),each=2)
plot(MA_rptr_df$R[order(MA_rptr_df$trait)]~ vect_x ,ylab="Repeatabillity",bty="n",las=1,xaxt="n",ylim=c(0,1.1),xlab="",type="n",yaxt="n")
arrows(vect_x,MA_rptr_df$lwr[order(MA_rptr_df$trait)],vect_x,MA_rptr_df$upr[order(MA_rptr_df$trait)],code=3,angle=90, length =.05)
points(MA_rptr_df$R[order(MA_rptr_df$trait)]~ vect_x,pch=21,bg= c("darkgreen","magenta")[as.numeric(as.factor(MA_rptr_df$pop))[order(MA_rptr_df$trait)]])
axis(side=1,at=c(4,12,20,28,36,44),labels=c("SF","SB","FS","FB","BS","BF"))
legend(2,1.2,c("MA lines - N2","MA lines - PB306"),bty="n",pch=c(16,16),col=c("darkgreen","magenta"),ncol=2)


save(list=ls(),file="~/PATH/TO/DIR/Cemee_Pop_WI/Repeatability.RData")
