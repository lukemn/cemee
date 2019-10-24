

rm(list=ls())


##### The model

nb_replicates=10
param_df <- data.frame(nb_gen=100,nb_gen_no_mut=rep(c(0),each=nb_replicates),nb_trait=6,nb_loci=100,K=1000,nb_off=10,lambda_pois=3,mu_rate_in=rep(c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7),each=nb_replicates),nb_loci_scaling_factor= 0,env_var=0,nb_init_alleles=16,max_alleles_per_loci= 100,nb_loci_scaling_factor_mut=NA,run_location="local",MA_line=NA,output_name=c("NA"), stringsAsFactors=FALSE, house_of_cards=NA)

param_list <- list()
for(i in 1:nrow(param_df)) param_list[[i]] <- param_df[i,]


library(parallel)
no_cores <- 30
clust <- makeCluster(no_cores)
output_list <- parLapply(clust, param_list , launch_main_model)
stopCluster(clust)

output_mat=NULL
for(i in 1:length(output_list)){
output_mat=rbind(output_mat,output_list[[i]])	
}
output_mat=output_mat+1
sav_output_mat=output_mat

output_mat=tapply(sav_output_mat[,1],param_df$mu,mean)
for(i in 2:ncol(sav_output_mat)) output_mat=cbind(output_mat,tapply(sav_output_mat[,i],param_df$mu,mean))
vect_mu=as.numeric(row.names(output_mat))
plot(output_mat[,1]~vect_mu,type="b",log="xy",ylim=c(1,5000),xlab=c("Mutation rate"),ylab="Number of new mutations",bty="n",xaxt="n")
axis(side=1,at=sort(c(1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)))
for(i in 2:5) points(output_mat[,i]~vect_mu,type="b",col=c("black","purple","blue","green","yellow","orange","red")[i])
legend(1e-7,5e+3,100*c(0.001,0.01,0.05,0.1,0.2,0.5,0.75)[1:5],lwd=2,col=c("black","purple","blue","green","yellow","orange","red")[1:5],title="Min. allele frequency (%)")

write.table(output_mat,file='~/PATH/TO/DIR/Simulations/output_Nb_mut.Rdata')

######################################################################################

launch_main_model <- function(arg){

library(matrixcalc)
library(MASS)
#library(matlib)
library(matrixStats)
library(MCMCglmm)
library(stats)
library(dae)
####
## Needed function
produce_Pop_gen_array <- function(nb_loci, K, nb_init_alleles = 16) {

	temp_freq_all <- rexp((nb_init_alleles) * nb_loci, 5) + 0.15
	temp_freq_all <- matrix(temp_freq_all[temp_freq_all < 0.9][1:I(nb_loci)], nb_loci, 
		((nb_init_alleles - 1)))

	both_chr_alleles <- 2 * round(K * temp_freq_all[, 1])

	temp_freq_all[, 1] <- round(2 * K * temp_freq_all[, 1])
	assigned_alleles <- 1

	while (assigned_alleles < c(nb_init_alleles - 1)) {
		assigned_alleles <- assigned_alleles + 1

		if (assigned_alleles == 2) 
			temp_freq_all[, 2] <- round(rowMaxs(cbind(rep(0, nb_loci), round(2 * K - 
				temp_freq_all[, 1]))) * temp_freq_all[, 2])
		if (assigned_alleles > 2) 
			temp_freq_all[, assigned_alleles] <- round(rowMaxs(cbind(rep(0, nb_loci), 
				round(2 * K - rowSums(temp_freq_all[, c(1:(assigned_alleles - 1))])))) * 
				temp_freq_all[, assigned_alleles])
	}
	temp_freq_all <- cbind(temp_freq_all, 2 * K - rowSums(temp_freq_all))


	## We have the frequency of all alleles, we will now create the CHR
	Pop_gen_array <- array(, c(K, nb_loci, 2))

	#for each loci
	for (i in 1:nb_loci) {
		## randomly place the alleles
		random_alleles_vector <- rep(c(1:nb_init_alleles), temp_freq_all[i, ])[sample(1:(2 * 
			K), 2 * K, replace = FALSE)]

		Pop_gen_array[, i, 1] <- random_alleles_vector[1:K]
		Pop_gen_array[, i, 2] <- random_alleles_vector[(1 + K):(2 * K)]
	}

	return(Pop_gen_array)
}



#Extract parameters from arg vector
nb_gen <- as.numeric(arg[1])
nb_gen_no_mut <- as.numeric(arg[2])
nb_trait <- as.numeric(arg[3])
nb_loci <- as.numeric(arg[4])
K <- as.numeric(arg[5])
nb_off <- as.numeric(arg[6])


mu_rate_in <- as.numeric(arg[8])
nb_init_alleles <- as.numeric(arg[11])
max_alleles_per_loci <- as.numeric(arg[12])

## Initialize path to save results

Males <- rep(FALSE, K)
Males[sample(c(1:K), I(K/2), replace = FALSE)] <- TRUE
mu <- mu_rate_in # mutation rate per locus

# Create the initial population

Pop_gen_array <- produce_Pop_gen_array(nb_loci, K, nb_init_alleles)

### Cycle on generations
for (gen_count in 1:nb_gen) {

	# All survived
	Survived <- rep(TRUE, K)

	# Mate the males with females
	nb_surviv_males <- sum(Males[Survived])
	mated_males <- c(1:K)[Males & Survived]
	# If not enough females, some males does not mate -> they die
	if (nb_surviv_males > I(sum(Survived)/2)) {
		nb_surviv_males = sum(!Males & Survived)
		mated_males <- sample(mated_males, nb_surviv_males)
	}

	mated_females <- sample(c(1:K)[!Males & Survived], nb_surviv_males, replace = FALSE)
	nb_produced_offs <- nb_off * nb_surviv_males

	#Produce the offsprings
	Offsprings <- array(, c(nb_produced_offs, nb_loci, 2))
	count_offs = 0

	# Offsprings from mated females	
	for (prod_off_mated in 1:nb_surviv_males) {
		for (vect_offs_prod in 1:nb_off) {
			count_offs <- count_offs + 1
			selected_chr1 <- round(runif(nb_loci,min=1,max=2))
			selected_chr2 <- round(runif(nb_loci,min=1,max=2))			
			Offsprings[count_offs, , ] <- cbind(diag(Pop_gen_array[mated_males[prod_off_mated], , selected_chr1]), 
				diag(Pop_gen_array[mated_females[prod_off_mated], , selected_chr2]))
		}
	}

	### Then replace the K parents by randomly sampled offsprings- 
	selected_offs <- sample(c(1:count_offs), K)
	Pop_gen_array <- Offsprings[selected_offs, , ]
	Males <- rep(FALSE, K)
	Males[sample(1:K, I(K/2), replace = FALSE)] <- TRUE # Sex of the individuals, TRUE if males

	if (gen_count > nb_gen_no_mut) {
		### New mutation
		mutated_alleles <- c(1:c(K * nb_loci * 2))[runif(K * nb_loci * 2) < mu]
		mutated_alleles <- mutated_alleles-1

		chr_mut <- floor(mutated_alleles/I(K * nb_loci)) + 1
		ind_mut <- floor(mutated_alleles%%I(K * nb_loci)/nb_loci) + 1
		loci_mut <- (mutated_alleles%%I(K * nb_loci))%%nb_loci +1
		new_allelic_values <- round(runif(length(mutated_alleles), min = (nb_init_alleles + 1), max = max_alleles_per_loci))
		
		for (vect_mutated_alleles in c(1:length(mutated_alleles))) {
			Pop_gen_array[ind_mut[vect_mutated_alleles], loci_mut[vect_mutated_alleles], chr_mut[vect_mutated_alleles]] <- new_allelic_values[vect_mutated_alleles]
		}
	}
	

} #End of the loop for the generations

over_Y_percent=function(x,Y=0.01){
	if(length(x)>16){
		return(sum(x[17:length(x)]>=Y*K*2))
	}else{
		return(0)
	}	
}

nb_mut_detected=NULL
for(i in c(0.001,0.01,0.05,0.1,0.2,0.5,0.75)){
nb_mut_detected=c(nb_mut_detected,sum(unlist(lapply(apply(Pop_gen_array,2,table), over_Y_percent,Y=i))))
}
return(nb_mut_detected)

}





