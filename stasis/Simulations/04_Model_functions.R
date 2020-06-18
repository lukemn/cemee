library(matrixcalc)
library(MASS)
#library(matlib)
library(matrixStats)
library(MCMCglmm)
library(stats)
library(dae)


angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}

Produce_Mut_Gen_arrays_for_G_simuls <- function(args) {

	source('~/Projets_Recherches/Celegans/G_matrix_manuscript/AmNat/Github_Rcodes_AmNat/Simulations/04_Model_functions.R', chdir = TRUE)
	A6140_mat <- read.table("~/Projets_Recherches/Celegans/G_matrix_manuscript/AmNat/Shared_files/G_mat_for_simulations/A6140.txt", sep = "\t")

	G_init_mat <- as.matrix(A6140_mat/2)

	nb_loci_scaling_factor <- as.numeric(args[1])
	lambda_pois <- as.numeric(args[2])
	nb_loci <- as.numeric(args[3])

	nb_trait <- 6 # Nb of traits
	K <- 150 # Number of lines

	nb_init_alleles <- 16
	max_alleles_per_loci <- 16


#### Create the population
Pop_gen_array <- produce_Pop_gen_array(nb_loci, K, nb_init_alleles)

#Associates mut effects
Mut_effect <- array(, c(max_alleles_per_loci, nb_loci, nb_trait))
Mut_effect <- assign_allelic_effects(max_alleles_per_loci, nb_loci, nb_trait, G_init_mat, 
		lambda_pois, nb_loci_scaling_factor)

	return(list(Pop_gen_array, Mut_effect, K))

}

##########
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


#####

assign_allelic_effects <- function(nb_alleles, nb_loci, nb_trait, G_Mat, lambda_pois, 
	scaling_mut_factor) {

	library(data.table)
	library(lcmix)

### We need to compute a relative effect size of alleles depending on the number
### of trait they affect
	temp <- rpois(n = 1e7, lambda = lambda_pois)
	temp=temp[temp>0 & temp<=nb_trait]
	temp = table(temp)
	temp=temp/sum(temp)

	rel_effect_size=(1/temp/(sum(1/temp))*nb_trait)


	Mut_effect <- array(, c(nb_alleles, nb_loci, nb_trait))
	tot_nb_alleles <- nb_alleles * nb_loci

	# generate the nb of traits affected by each alleles
	nb_alleles_affected <- 0
	vect_nb_mutated_traits <- NULL
	while (nb_alleles_affected < tot_nb_alleles) {
		nb_traits_temp <- rpois(n = tot_nb_alleles, lambda = lambda_pois)
		vect_nb_mutated_traits <- c(vect_nb_mutated_traits, nb_traits_temp[nb_traits_temp > 
			0 & nb_traits_temp <= nb_trait])
		nb_alleles_affected = length(vect_nb_mutated_traits)
	}
	vect_nb_mutated_traits <- table(factor(vect_nb_mutated_traits[1:tot_nb_alleles], 
		levels = 1:nb_trait))

	### mutations that affects only one trait
	mat_allelic_effects <- matrix(rep(0, nb_trait * tot_nb_alleles), ncol = nb_trait)
	vect_sampled_traits <- sample(1:nb_trait, vect_nb_mutated_traits[[1]], replace = TRUE)

	vect_allelic_effects <- rnorm(vect_nb_mutated_traits[1], mean = rep(0, 6), diag(G_Mat)[vect_sampled_traits]) * rel_effect_size[1]

	mat_allelic_effects[1:vect_nb_mutated_traits[1], vect_sampled_traits] <- vect_allelic_effects
	k_count = vect_nb_mutated_traits[1]

	### mutations that affects two or more traits
	for (i_nb_traits in 2:nb_trait) {

		#vect_nb_mutated_traits[i_nb_traits]
		vect_sampled_traits <- NULL
		if(vect_nb_mutated_traits[[i_nb_traits]]>0){
		for (i in 1:vect_nb_mutated_traits[[i_nb_traits]]) vect_sampled_traits <- c(vect_sampled_traits, 
			paste(sort(sample(1:nb_trait, i_nb_traits, replace = FALSE)), collapse = " "))

		vect_sampled_traits <- table(vect_sampled_traits)
		list_T <- list()
		for (j_nb_traits in 1:i_nb_traits) list_T[[j_nb_traits]] <- as.numeric(tstrsplit(names(vect_sampled_traits), 
			" ", fixed = TRUE)[[j_nb_traits]])

		for(j_nb_combin in 1:length(list_T[[1]])) {

			temp_vect_traits = NULL
			for (j_temp_traits in 1:length(list_T)) temp_vect_traits <- c(temp_vect_traits, 
				list_T[[j_temp_traits]][j_nb_combin])

			temp_allelic_effects <- rmvnorm(n = vect_sampled_traits[j_nb_combin], mean = rep(0, 
				i_nb_traits), cov = matrix(G_Mat[temp_vect_traits, temp_vect_traits], 
				ncol = i_nb_traits)) * rel_effect_size[i_nb_traits]

			mat_allelic_effects[I(k_count + 1):I(k_count + vect_sampled_traits[j_nb_combin]), 
				temp_vect_traits] <- temp_allelic_effects
			k_count = k_count + vect_sampled_traits[j_nb_combin]

		}
	}
}
	#Randomize the order of the allelic effects
	mat_allelic_effects <- mat_allelic_effects[sample(1:tot_nb_alleles, tot_nb_alleles, 
		replace = FALSE), ]
	#Assign alleles for each loci
	k_count = 0
	for (i_affect_all in 1:nb_alleles) {
		Mut_effect[i_affect_all, , ] = mat_allelic_effects[(k_count + 1):(k_count + 
			nb_loci), ]
		k_count = k_count + nb_loci
	}
	Mut_effect = Mut_effect * scaling_mut_factor
	return(Mut_effect)
}

#####

compute_G_distances <- function(G_mat_ref, G_exp_mat) {

	min_i = 1
	max_i = 6
	vect_angles <- NULL
	vect_ratios <- NULL
	vect_summed_eigVal <- NULL

	for (i1 in c(min_i:max_i)) {
		for (i2 in c(min(max_i, (i1 + 1)):max_i)) {
			if (i1 != i2) {
				A <- G_exp_mat[c(i1, i2), c(i1, i2)]
				B <- G_mat_ref[c(i1, i2), c(i1, i2)]

				eigVal_A <- eigen(A)$values
				eigVec_A <- eigen(A)$vectors
				eigVal_B <- eigen(B)$values
				eigVec_B <- eigen(B)$vectors

				vect_angles <- c(vect_angles, angle_eigenV(eigVec_B[1, ], eigVec_A[1, 
					]))
				vect_angles[vect_angles > I(pi/2)] <- pi - vect_angles[vect_angles > I(pi/2)]
				vect_ratios <- c(vect_ratios, eigVal_A[1]/eigVal_A[2] - eigVal_B[1]/eigVal_B[2])
				vect_summed_eigVal <- c(vect_summed_eigVal, (sum(eigVal_A) - sum(eigVal_B))/sum(eigVal_B))

			}
		}
	}

	return(list(sum(vect_angles), sum(vect_ratios), sum(vect_summed_eigVal)))

}


##########





process_output_from_simulations <- function(list_for_processing) {


	output_item <- list_for_processing[[1]]

	source('~/Projets_Recherches/Celegans/G_matrix_manuscript/AmNat/Github_Rcodes_AmNat/Simulations/04_Model_functions.R', chdir = TRUE)


	#reformat the arrays
	list_gen_mut_arrays <- list()
	for (j in 1:2) list_gen_mut_arrays[[j]] <- list(output_item$Pop_gen_array_sav[[j]], 
		output_item$Mut_effect_sav[[j]], 150)

	temp_list <- lapply(list_gen_mut_arrays, produce_G_from_phen_array)
	return(temp_list)
}




assign_allelic_effects <- function(nb_alleles, nb_loci, nb_trait, G_Mat, lambda_pois, 
	scaling_mut_factor) {

	library(data.table)
	library(lcmix)

### We need to compute a relative effect size of alleles depending on the number
### of trait they affect
	temp <- rpois(n = 1e7, lambda = lambda_pois)
	temp=temp[temp>0 & temp<=nb_trait]
	temp = table(temp)
	temp=temp/sum(temp)

	rel_effect_size=(1/temp/(sum(1/temp))*nb_trait)


	Mut_effect <- array(, c(nb_alleles, nb_loci, nb_trait))
	tot_nb_alleles <- nb_alleles * nb_loci

	# generate the nb of traits affected by each alleles
	nb_alleles_affected <- 0
	vect_nb_mutated_traits <- NULL
	while (nb_alleles_affected < tot_nb_alleles) {
		nb_traits_temp <- rpois(n = tot_nb_alleles, lambda = lambda_pois)
		vect_nb_mutated_traits <- c(vect_nb_mutated_traits, nb_traits_temp[nb_traits_temp > 
			0 & nb_traits_temp <= nb_trait])
		nb_alleles_affected = length(vect_nb_mutated_traits)
	}
	vect_nb_mutated_traits <- table(factor(vect_nb_mutated_traits[1:tot_nb_alleles], 
		levels = 1:nb_trait))

	### mutations that affects only one trait
	mat_allelic_effects <- matrix(rep(0, nb_trait * tot_nb_alleles), ncol = nb_trait)
	vect_sampled_traits <- sample(1:nb_trait, vect_nb_mutated_traits[[1]], replace = TRUE)

	vect_allelic_effects <- rnorm(vect_nb_mutated_traits[1], mean = rep(0, 6), diag(G_Mat)[vect_sampled_traits]) * rel_effect_size[1]

	mat_allelic_effects[1:vect_nb_mutated_traits[1], vect_sampled_traits] <- vect_allelic_effects
	k_count = vect_nb_mutated_traits[1]

	### mutations that affects two or more traits
	for (i_nb_traits in 2:nb_trait) {

		#vect_nb_mutated_traits[i_nb_traits]
		vect_sampled_traits <- NULL
		if(vect_nb_mutated_traits[[i_nb_traits]]>0){
		for (i in 1:vect_nb_mutated_traits[[i_nb_traits]]) vect_sampled_traits <- c(vect_sampled_traits, 
			paste(sort(sample(1:nb_trait, i_nb_traits, replace = FALSE)), collapse = " "))

		vect_sampled_traits <- table(vect_sampled_traits)
		list_T <- list()
		for (j_nb_traits in 1:i_nb_traits) list_T[[j_nb_traits]] <- as.numeric(tstrsplit(names(vect_sampled_traits), 
			" ", fixed = TRUE)[[j_nb_traits]])

		for(j_nb_combin in 1:length(list_T[[1]])) {

			temp_vect_traits = NULL
			for (j_temp_traits in 1:length(list_T)) temp_vect_traits <- c(temp_vect_traits, 
				list_T[[j_temp_traits]][j_nb_combin])

			temp_allelic_effects <- rmvnorm(n = vect_sampled_traits[j_nb_combin], mean = rep(0, 
				i_nb_traits), cov = matrix(G_Mat[temp_vect_traits, temp_vect_traits], 
				ncol = i_nb_traits)) * rel_effect_size[i_nb_traits]

			mat_allelic_effects[I(k_count + 1):I(k_count + vect_sampled_traits[j_nb_combin]), 
				temp_vect_traits] <- temp_allelic_effects
			k_count = k_count + vect_sampled_traits[j_nb_combin]

		}
	}
}
	#Randomize the order of the allelic effects
	mat_allelic_effects <- mat_allelic_effects[sample(1:tot_nb_alleles, tot_nb_alleles, 
		replace = FALSE), ]
	#Assign alleles for each loci
	k_count = 0
	for (i_affect_all in 1:nb_alleles) {
		Mut_effect[i_affect_all, , ] = mat_allelic_effects[(k_count + 1):(k_count + 
			nb_loci), ]
		k_count = k_count + nb_loci
	}
	Mut_effect = Mut_effect * scaling_mut_factor
	return(Mut_effect)
}


produce_G_from_phen_array <- function(List_pop_mut_array, K = 100, nb_rep_lines = 2, 
	env_var = 0) {
	library(MCMCglmm)
	library(data.table)
	Pop_gen_array <- List_pop_mut_array[[1]]
	Mut_effect <- List_pop_mut_array[[2]]
	K <- min(List_pop_mut_array[[3]], dim(Pop_gen_array)[1])

	nb_loci <- dim(Pop_gen_array)[2]
	nb_trait = dim(Mut_effect)[3]

	## Estimate lambda pois
	#k_count = 0
	#summed_nb_eff = 0
	#for (i in 1:dim(Mut_effect)[1]) {
	#	summed_nb_eff = summed_nb_eff + sum(Mut_effect[i, , ] != 0)
	#	k_count = k_count + nrow(Mut_effect[i, , ])
	#}
	#lambda_pois_est <- summed_nb_eff/k_count

	## Phenotypic values of individuals
	Pop_phen_array <- NULL
	for (i in sample(1:nrow(Pop_gen_array),K)){
	rdm_chr=sample(1:2,1)
		temp_df = data.frame(Line = rep(paste0("L", i), nb_rep_lines))

		for (vect_t in c(1:nb_trait)) {
			temp_df$temp = sum(diag(Mut_effect[Pop_gen_array[i, , rdm_chr], , vect_t]) + diag(Mut_effect[Pop_gen_array[i,,rdm_chr], , vect_t]))
			names(temp_df)[ncol(temp_df)] <- paste0("T", vect_t)
		}
		Pop_phen_array = rbind(Pop_phen_array, temp_df)
	}
	#for (i_adjust_P in 1:nrow(Pop_phen_array)) Pop_phen_array[i_adjust_P, 2:7] <- Pop_phen_array[i_adjust_P,2:7]

	############################ Estimate the G Matrix
	vect_P_traits = paste0("T", 1:6)
	phen.var = diag(nb_trait) * diag(var(subset(Pop_phen_array, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/2, n = nb_trait)), R = list(V = phen.var/2, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T1, T2, T3, T4, T5, T6)) ~ trait - 1, random = ~us(trait):Line, 
		rcov = ~us(trait):units, family = rep("gaussian", nb_trait), data = Pop_phen_array, 
		prior = prior_mod, verbose = FALSE)
	post_dist = posterior.mode(model_MCMC$VCV)

	return(list(G1_mat = matrix(post_dist[1:nb_trait^2], nb_trait, nb_trait), R_mat = matrix(post_dist[(nb_trait^2 + 
		1):(2 * nb_trait^2)], nb_trait, nb_trait), VCV_mat = model_MCMC$VCV))
}






produce_Phenotypes_from_phen_array <- function(List_pop_mut_array, K = 100, nb_rep_lines = 2, 
	env_var = 0) {

	Pop_gen_array <- List_pop_mut_array[[1]]
	Mut_effect <- List_pop_mut_array[[2]]
	nb_loci <- dim(Pop_gen_array)[2]
	nb_trait = dim(Mut_effect)[3]
	## Phenotypic values of individuals
	Pop_phen_array <- NULL
	for (i in c(1:K)) {

		temp_df = data.frame(Line = rep(paste0("L", i), nb_rep_lines))

		for (vect_t in c(1:nb_trait)) {
			temp_df$temp = sum(diag(Mut_effect[Pop_gen_array[i, , 1], , vect_t]) + diag(Mut_effect[Pop_gen_array[i, 
				, 2], , vect_t])) + rnorm(nb_rep_lines, 0, sqrt(env_var))
			names(temp_df)[ncol(temp_df)] <- paste0("T", vect_t)
		}

		Pop_phen_array = rbind(Pop_phen_array, temp_df)

	}
	return(Pop_phen_array)
}



process_output_from_simulations <- function(list_for_processing) {


	output_item <- list_for_processing[[1]]

	source('~/Projets_Recherches/Celegans/G_matrix_manuscript/AmNat/Github_Rcodes_AmNat/Simulations/Model_functions.R', chdir = TRUE)

	#reformat the arrays
	list_gen_mut_arrays <- list()
	for (j in 1:2) list_gen_mut_arrays[[j]] <- list(output_item$Pop_gen_array_sav[[j]], 
		output_item$Mut_effect_sav[[j]], 150)

	temp_list <- lapply(list_gen_mut_arrays, produce_G_from_phen_array)
	return(temp_list)
}


produce_G_from_phen_array_v2 <- function(List_pop_mut_array, nb_rep_lines = 2) {

	Pop_gen_array <- List_pop_mut_array[[1]]
	Mut_effect <- List_pop_mut_array[[2]]
	K <- min(List_pop_mut_array[[5]], dim(Pop_gen_array)[1])

	nb_loci <- dim(Pop_gen_array)[2]
	nb_trait = dim(Mut_effect)[3]

	## Estimate lambda pois
	k_count = 0
	summed_nb_eff = 0
	for (i in 1:dim(Mut_effect)[1]) {
		summed_nb_eff = summed_nb_eff + sum(Mut_effect[i, , ] != 0)
		k_count = k_count + nrow(Mut_effect[i, , ])
	}
	lambda_pois_est <- summed_nb_eff/k_count


	## Phenotypic values of individuals
	Pop_phen_array <- NULL
	for (i in c(1:K)) {

		temp_df = data.frame(Line = rep(paste0("L", i), nb_rep_lines))
		for (vect_t in c(1:nb_trait)) {
			#Inbred lines
temp_df$temp = sum(diag(Mut_effect[Pop_gen_array[i, , 1], , vect_t]) + diag(Mut_effect[Pop_gen_array[i, 
				, 1], , vect_t]))


			names(temp_df)[ncol(temp_df)] <- paste0("T", vect_t)
		}

		Pop_phen_array = rbind(Pop_phen_array, temp_df)

	}

	############################ Estimate the G Matrix
	vect_P_traits = paste0("T", 1:6)
;	phen.var = diag(nb_trait) * diag(var(subset(Pop_phen_array, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, 
		n = nb_trait)), R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T1, T2, T3, T4, T5, T6)) ~ trait - 1, random = ~us(trait):Line + 
		us(trait):block, rcov = ~us(trait):units, family = rep("gaussian", nb_trait), 
		data = Pop_phen_array, prior = prior_mod, verbose = FALSE)
	post_dist = posterior.mode(model_MCMC$VCV)

	return(list(G1_mat = matrix(post_dist[1:nb_trait^2], nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 
		1):(2 * nb_trait^2)], nb_trait, nb_trait), R_mat = matrix(post_dist[(2 * nb_trait^2 + 
		1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV))

}



compute_G_distances <- function(G_mat_ref, G_exp_mat) {

	min_i = 1
	max_i = 6
	vect_angles <- NULL
	vect_ratios <- NULL
	vect_summed_eigVal <- NULL

	for (i1 in c(min_i:max_i)) {
		for (i2 in c(min(max_i, (i1 + 1)):max_i)) {
			if (i1 != i2) {
				A <- G_exp_mat[c(i1, i2), c(i1, i2)]
				B <- G_mat_ref[c(i1, i2), c(i1, i2)]

				eigVal_A <- eigen(A)$values
				eigVec_A <- eigen(A)$vectors
				eigVal_B <- eigen(B)$values
				eigVec_B <- eigen(B)$vectors

				vect_angles <- c(vect_angles, angle_eigenV(eigVec_B[1, ], eigVec_A[1, 
					]))
				vect_angles[vect_angles > I(pi/2)] <- pi - vect_angles[vect_angles > I(pi/2)]
				vect_ratios <- c(vect_ratios, eigVal_A[1]/eigVal_A[2] - eigVal_B[1]/eigVal_B[2])
				vect_summed_eigVal <- c(vect_summed_eigVal, (sum(eigVal_A) - sum(eigVal_B))/sum(eigVal_B))

			}
		}
	}

	return(list(sum(vect_angles), sum(vect_ratios), sum(vect_summed_eigVal)))

}

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}

