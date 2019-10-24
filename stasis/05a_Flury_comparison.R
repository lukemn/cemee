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

load("~/PATH/TO/DIR/Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData")

#Export tables for the Flury comparisons
data_for_Flury= data_for_G_estimation[order(data_for_G_estimation$population, data_for_G_estimation$pop_label),]

for(i in paste0(vect_P_traits,"_pred_02")){
data_for_Flury[,i]= (data_for_Flury[,i]-mean(data_for_Flury[,i]))
}

vect_p=unique(data_for_Flury$population)[c(1,3,5,7,2,4,6)]

for(i in 1:6){
	for(j in (i+1):7){

write.table(subset(data_for_Flury,population%in% vect_p[c(i,j)])[,c("population","pop_label",paste0(vect_P_traits,"_pred_02"))],file=paste0("~/PATH/TO/DIR/Matrix_comparison/tables_original/all_traits_",vect_p[i],"_",vect_p[j],".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
		}
	}

## Then I wrote an "expect" script that can be ran from the terminal in MAC OS
## cd ~/PATH/TO/DIR/Matrix_comparison
## expect ./expect_script_original.sh


### We can extract the p values from all tests

setwd("~/PATH/TO/DIR/Matrix_comparison/results_original")

all_results_files = list.files()
length(all_results_files) # 21

splitted_char = tstrsplit(all_results_files, "_")
splitted_char2 = tstrsplit(splitted_char[[4]], ".txt")[[1]]

flury_results = data.frame(pop1= splitted_char[[3]],pop2= splitted_char2,G1=NA,G2=NA,G3=NA,G4=NA,G5=NA,G6=NA,G7=NA)
rm(splitted_char, splitted_char2)

for(i in 1:length(all_results_files)){
temp= readLines(all_results_files[i])

add_index =ifelse(substring(temp[54],1,7)=="Bending",1,0)
flury_results$G1[i] =as.numeric(tstrsplit(temp[62-add_index],"p = ")[[2]])
flury_results$G2[i] =as.numeric(tstrsplit(temp[67-add_index],"p = ")[[2]])
flury_results$G3[i] =as.numeric(tstrsplit(temp[72-add_index],"p = ")[[2]])

flury_results$G4[i] =as.numeric(tstrsplit(temp[77-add_index],"p = ")[[2]])
flury_results$G5[i] =as.numeric(tstrsplit(temp[82-add_index],"p = ")[[2]])
flury_results$G6[i] =as.numeric(tstrsplit(temp[87-add_index],"p = ")[[2]])
flury_results$G7[i] =as.numeric(tstrsplit(temp[92-add_index],"p = ")[[2]])
}


thresh = .05
flury_results$final_label = NA

flury_results$final_label[rowSums(flury_results[,3:9]< thresh)==7]= "diff"

flury_results$final_label[flury_results[,9]> thresh]= "CPC1"
flury_results$final_label[flury_results[,8]> thresh]= "CPC2"
flury_results$final_label[flury_results[,7]> thresh]= "CPC3"
flury_results$final_label[flury_results[,6]> thresh]= "CPC4"
flury_results$final_label[flury_results[,5]> thresh]= "CPC"
flury_results$final_label[flury_results[,4]> thresh]= "prop"
flury_results$final_label[flury_results[,3]> thresh]= "equal"

subset(flury_results,is.na(final_label))

v_CA50 = c("CA150","CA250","CA350")
v_CA100 = c("CA1100","CA2100","CA3100")
flury_results$comp_type = "CA50 - CA100"
flury_results[flury_results$pop1=="A6140" & flury_results$pop2%in% v_CA50,]$comp_type="A6 - CA50"
flury_results[flury_results$pop1=="A6140" & flury_results$pop2%in% v_CA100,]$comp_type="A6 - CA100"

flury_results[flury_results$pop1%in%v_CA50 & flury_results$pop2%in% v_CA50,]$comp_type="CA50 - CA50"
flury_results[flury_results$pop1%in%v_CA100 & flury_results$pop2%in% v_CA100,]$comp_type="CA100 - CA100"

# A table can then be saved
write.table(flury_results,file="~/PATH/TO/DIR/Matrix_comparison/Flury_table.txt",quote=FALSE,row.names=FALSE,sep="\t")

