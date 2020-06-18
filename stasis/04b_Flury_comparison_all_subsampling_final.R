rm(list = ls())
gc()
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

### Compute a null distribution of LR for the Flury comparison (Table S2)

load('~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData')

#Export tables for the Flury comparisons
data_for_Flury= data_for_G_estimation[order(data_for_G_estimation$population, data_for_G_estimation$pop_label),]

for(i in paste0(vect_P_traits,"_pred_02")){
data_for_Flury[,i]= (data_for_Flury[,i]-mean(data_for_Flury[,i]))
}
data_for_Flury$population=as.factor(as.character(data_for_Flury$population))

output_final_list=list()
output_final_df=NULL
sign=NULL

# We will build a distribution of statistics/p-values for the comparisons
# -> one for the comparision A6140/CA[1-3] and the other for the comparison
# CA[1-3]/CA[1-3]

#Number of lines in the A6140
nb_popA6=length(unique(subset(data_for_Flury,population=="A6140")$pop_label))

nb_popCA=mean(table(unique(data_for_Flury[,c("population","pop_label")])$population)[2:7])

#If we sample from all

unique_lines_label=unique(data_for_Flury$pop_label)
unique_lines_label_noA6=unique(subset(data_for_Flury,population!="A6140")$pop_label)
unique_lines_label_A6=unique(subset(data_for_Flury,population=="A6140")$pop_label)

length(unique_lines_label_noA6) # 327
length(unique_lines_label) 	# 516

flury_results_dist=list()

for(vect_k_comp in 1:2){
#vect_k_comp=1
flury_results_temp_LR = data.frame(G1=NA,G2=NA,G3=NA,G4=NA,G5=NA,G6=NA,G7=NA)

for(i in 1:1000){
	if(i%%10==0) print(i)
	if(vect_k_comp==1) unique_lines = unique_lines_label
	if(vect_k_comp==2) unique_lines = unique_lines_label_noA6
	#unique_lines = unique_lines#_label_A6	

	set1 = sample(unique_lines, nb_popA6,replace=TRUE)
	set2 = sample(unique_lines, nb_popCA,replace=TRUE)
	set3 = sample(unique_lines, nb_popCA,replace=TRUE)	
	
	d1=subset(data_for_Flury,pop_label%in%set1);d1$population="set1"
	d2=subset(data_for_Flury,pop_label%in%set2);d2$population="set2"
	d3=subset(data_for_Flury,pop_label%in%set3);d3$population="set3"
	
	data_for_Flury_temp1=rbind(d1,d2)
	data_for_Flury_temp2=rbind(d2,d3)

data_for_Flury_temp1 = data_for_Flury_temp1[order(data_for_Flury_temp1$population, data_for_Flury_temp1$pop_label),]
data_for_Flury_temp2 = data_for_Flury_temp2[order(data_for_Flury_temp2$population, data_for_Flury_temp2$pop_label),]

if(vect_k_comp==1) data_for_Flury_temp= data_for_Flury_temp1
if(vect_k_comp==2) data_for_Flury_temp= data_for_Flury_temp2
	write.table(data_for_Flury_temp[,c("population","pop_label",paste0(vect_P_traits,"_pred_02"))],file="~/PATH/TO/DIR/Flury_tables_all/all_traits_rdm.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

system("rm ~/PATH/TO/DIR/Flury_results_all/*")

output_script=system("~/PATH/TO/DIR/Flury_scripts/expect_script_original_all.sh > temp_to_delete.txt")

all_output_files=list.files("~/PATH/TO/DIR/Flury_results_all/")

### We can extract the test distribution values from the file

temp=readLines(paste0("~/PATH/TO/DIR/Flury_results_all/",all_output_files))
add_index =ifelse(substring(temp[54],1,7)=="Bending",1,0)

flury_results_temp_LR=rbind(flury_results_temp_LR,c(
as.numeric(tstrsplit(temp[60-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[65-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[70-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[75-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[80-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[85-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[90-add_index],"estimates: ")[[2]])))
}

flury_results_temp_LR =subset(flury_results_temp_LR,!is.na(G1))
flury_results_dist[[vect_k_comp]]= flury_results_temp_LR
}

## The true values

list_true_results=list.files("~/PATH/TO/DIR/Matrix_comparison/results_original/")

list_true_results_pop1=tstrsplit(list_true_results,"_")[[3]]
list_true_results_pop2=tstrsplit(tstrsplit(list_true_results,"_")[[4]],".txt")[[1]]

# Results in order
pop_cp1=rep(c("A6140","CA150","CA250","CA350","CA150","CA250","CA1100","CA2100"),c(6,1,1,1,2,1,2,1))
pop_cp2=c("CA150","CA250","CA350","CA1100","CA2100","CA3100","CA1100","CA2100","CA3100","CA250","CA350","CA350","CA2100","CA3100","CA3100")

flury_results_true=data.frame(G1=NA,G2=NA,G3=NA,G4=NA,G5=NA,G6=NA,G7=NA)
for(k in 1:length(pop_cp1)){

temp=readLines(paste0("~/PATH/TO/DIR/Matrix_comparison/results_original/",list_true_results[list_true_results_pop1==pop_cp1[k] & list_true_results_pop2==pop_cp2[k]]))
add_index =ifelse(substring(temp[54],1,7)=="Bending",1,0)

flury_results_true=rbind(flury_results_true,c(
as.numeric(tstrsplit(temp[60-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[65-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[70-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[75-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[80-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[85-add_index],"estimates: ")[[2]]),
as.numeric(tstrsplit(temp[90-add_index],"estimates: ")[[2]])))

}
flury_results_true=subset(flury_results_true,!is.na(flury_results_true$G1))

flury_results_true[1:6,]
sign_levels=apply(flury_results_dist[[1]],2,function(x){ sort(x)[(.95)*length(x)]})

sign_matrix=matrix(NA,15,7)

for(i in 1:7){
	for(j in 1:15){
		k=ifelse(j<=6,1,2)
		sign_matrix[j,i]		 =
		1-(sum(flury_results_true[j,i] > flury_results_dist[[k]][,i]) / nrow(flury_results_dist[[1]]))
		
	}}

write.table(sign_matrix,file='~/PATH/TO/DIR/Flury_S2.txt',quote=FALSE,sep="\t")

