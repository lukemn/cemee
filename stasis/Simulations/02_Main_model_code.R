rm(list=ls())

source('Simulations/Model_functions.R', chdir = TRUE)

nb_replicates <- 20
final_estimate = .1086921

## The order of the columns is critical!
param_df <- data.frame(nb_gen=100,nb_gen_no_mut=rep(c(101),each=nb_replicates),nb_trait=6,nb_loci=100,K=1000,nb_off=10,lambda_pois=3,mu_rate_in=0,nb_loci_scaling_factor= final_estimate,env_var=0,nb_init_alleles=16,max_alleles_per_loci= 16,nb_loci_scaling_factor_mut=NA,run_location="local",MA_line=NA,output_name=c("NA"), stringsAsFactors=FALSE)

param_list <- list()
for(i in 1:nrow(param_df)) param_list[[i]] <- param_df[i,]

##
library(parallel)
no_cores <- 3#detectCores()-2
clust <- makeCluster(no_cores)
output_list <- parLapply(clust, param_list , launch_main_model)
stopCluster(clust)


G_mat_start =list()
G_mat_final=list()
Gmat_list =list()
model_VCV_init=list()
model_VCV_evolved=list()


list_for_output_processing <- list()
for(i in 1:length(output_list)){
	 list_for_output_processing[[i]] <- list(output_list[[i]])
	 }
no_cores=3
clust <- makeCluster(no_cores)
list_model_output <- parLapply(clust,list_for_output_processing,process_output_from_simulations)
stopCluster(clust)


G_dist_all=NULL
for(kkk in 1:length(list_model_output)){

model_VCV_init[[kkk]] <- list_model_output[[kkk]][[1]]
model_VCV_evolved[[kkk]] <- list_model_output[[kkk]][[2]]

Gmat_list[[kkk]] <- list( model_VCV_init[[kkk]]$G1_mat/2 , model_VCV_evolved[[kkk]]$G1_mat/2)
G_dist <- data.frame(compute_G_distances(Gmat_list[[kkk]][[1]], Gmat_list[[kkk]][[2]]))
names(G_dist) <- c("angle","ratio","DiffeigenV")

G_dist$mut <- output_list[[kkk]]$nb_gen_no_mut>0
G_dist <- cbind(output_list[[kkk]]$arg,G_dist)
G_dist_all=rbind(G_dist_all, G_dist)

G_mat_start[[kkk]]=Gmat_list[[kkk]][[1]]
G_mat_final[[kkk]]=Gmat_list[[kkk]][[2]]

}

G_init_mat_path="~/PATH/TO/DIR/Cemee_Pop_WI/G_mat_for_simulations/A6140.txt"
A6140_mat <-read.table(G_init_mat_path)/2
E1_A6=eigen(A6140_mat)$vectors[,1]

G_dist_all$main_Angle <- NA

for(kkk in 1:nrow(G_dist_all)){
e1 <- eigen(Gmat_list[[kkk]][[1]])$vectors[,1]
e2 <- eigen(Gmat_list[[kkk]][[2]])$vectors[,1]

G_dist_all$main_Angle[kkk] <- angle_eigenV(e1,e2)

G_dist_all$main_Start[kkk] <- angle_eigenV(E1_A6,e1)
G_dist_all$main_End[kkk] <- angle_eigenV(E1_A6,e2)

}

save(list=ls(),file="~/PATH/TO/DIR/Simulations/Main_model_output_list_with_metadata.RData")


####################################################################################
####################################################################################

launch_main_model <- function(arg){

#Extract parameters from arg vector
nb_gen <- as.numeric(arg[1])
nb_gen_no_mut <- as.numeric(arg[2])
nb_trait <- as.numeric(arg[3])
nb_loci <- as.numeric(arg[4])
K <- as.numeric(arg[5])
nb_off <- as.numeric(arg[6])
lambda_pois <- as.numeric(arg[7])
mu_rate_in <- as.numeric(arg[8])
nb_loci_scaling_factor <- as.numeric(arg[9])
env_var <- as.numeric(arg[10])
nb_init_alleles <- as.numeric(arg[11])
max_alleles_per_loci <- as.numeric(arg[12])
nb_loci_scaling_factor_mut <- as.numeric(arg[13])
#run_location <- as.character(arg[14])
MA_line <-  as.character(arg[15])
output_name <-  as.character(arg[16])


path_to_functions='Simulations/Model_functions.R'
G_init_mat_path="~/PATH/TO/DIR/Cemee_Pop_WI/G_mat_for_simulations/A6140.txt"
A6140_mat <-read.table(G_init_mat_path)

source(path_to_functions, chdir = TRUE)

## Initialize path to save results
Mut_effect_sav <- list()
Pop_gen_array_sav <- list()
vect_sav <- 0

Males <- rep(FALSE, K)
Males[sample(c(1:K), I(K/2), replace = FALSE)] <- TRUE

G_init_mat <- as.matrix(A6140_mat/2)
mu <- mu_rate_in # mutation rate per locus

scaling_mut_G_factor <-  nb_loci_scaling_factor
theta <- rep(0, nb_trait)

# Create the initial population

Pop_gen_array <- produce_Pop_gen_array(nb_loci, K, nb_init_alleles)


#Assign allele effect on traits for each loci:alleles
		Mut_effect <- array(, c(max_alleles_per_loci, nb_loci, nb_trait)) # max_alleles_per_loci, each affecting all nb_trait		
		Mut_effect[1: nb_init_alleles,,] <- assign_allelic_effects(nb_init_alleles,nb_loci, nb_trait , G_init_mat ,lambda_pois, scaling_mut_G_factor)

#### Save the initial state

vect_sav= vect_sav+1
Mut_effect_sav [[vect_sav]]   <- Mut_effect
Pop_gen_array_sav[[vect_sav]] <- Pop_gen_array


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

vect_sav= vect_sav+1
Mut_effect_sav [[vect_sav]]   <- Mut_effect
Pop_gen_array_sav[[vect_sav]] <- Pop_gen_array


output_list <- list(arg ,Mut_effect_sav, Pop_gen_array_sav, nb_gen,nb_gen_no_mut,nb_trait,nb_loci,
K,nb_off,G_init_mat_path,lambda_pois,mu_rate_in,nb_loci_scaling_factor,
env_var,nb_init_alleles,max_alleles_per_loci)
names(output_list) <- c("arg" ,"Mut_effect_sav", "Pop_gen_array_sav", "nb_gen","nb_gen_no_mut","nb_trait","nb_loci",
"K","nb_off","G_init_mat_path","lambda_pois","mu_rate_in","nb_loci_scaling_factor",
"env_var","nb_init_alleles","max_alleles_per_loci")

return(output_list)

}




