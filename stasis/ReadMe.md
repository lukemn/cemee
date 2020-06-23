Supplementary code for [Selection and drift determine phenotypic stasis despite genetic divergence](https://www.biorxiv.org/content/10.1101/778282v1.full).

**01_Analysis_MA_Lines.R**
Loads the transition rates for the mutation accumulation lines (Final_merged_data_MA_Lines.txt), performs data normalization and estimates the **M**-matrices for each ancestor (N2, PB306), then outputs the matrices in different file formats.

**02_Analysis_Cemee_Pop_WI.R**
Loads the transition rates for the CeMEE lines (Final_merged_data_NGM.txt) , performs data normalization and estimates the **M**-matrices for each of the ancestors (N2, PB306), then outputs the matrices in different file formats.

**02a_State_frequency.R**
Loads the transition rates (Final_merged_data_NGM.txt) and the average observed frequencies of each state (State_freq_NGM.txt) and compares them, producing Figure B3.

**02b_Beta_distributions_with_freq.R**
Pre-processes the population data to produce the analysis and Figures in files [02c-j] (hermaphrodites).

**02c_Beta_distributions_with_freq_males.R**
Pre-processes the population data (males).

**02d_Figure2A.R**
Produces Figure 2A.
	
**02e_Repeatability.R**
Produces repeatability analysis and Figures 2D and B5.

**02f_Figure2BC**
Produces Figures 2B and C.

**02g_Male_frequency_and_Inbreeding**
Analyzes the effect of male frequency in the outbred population and the effect of inbreeding depression on transition rates. Produces Figures 3A and 3B.

**02h_Figure_4**
Produces Figure 4.

**02i_Figure5**
Produces Figure 5.

**03a_Gmatrix_WI.R**
Perform data normalization and estimate the **G**-matrix of the founders.

**03b_Gamma_matrix_estimation_with_betas_FigB6**
Perform estimation of the linear selection coefficients &beta; together with the $\gamma$-matrix from the fertility data. Produce Figure B6.

**03c_Gamma_matrix_estimation.R**
Estimate the $\gamma$-matrix from the fertility data. Compute the null distribution of selection on the rotated trait $y_i$ using mvnpermute and the *A* matrices of the CeMEE lines. Produce Figures 7A and B7.

**03d_Projection_canonical_axis.R**
Perform the canonical analysis of the transition rates based on the $\gamma$-matrix.

**03e_Figure_7B.R**
Produce Figure 7B from the canonical analysis.

**03f_Figure_B8.R**
Produce Figure B8 from the **M**-Matrices and the canonical axis.

**04a_Flury_comparison.R**
Process the Flury comparison of the **G**-matrices from the raw transition rates data using the “cpcrandMAC” executable (called from *R* using an expect script). Save the results for Table C8.

**04b_Flury_comparison_all_subsampling_final.R**
Produce the null distribution of LR statistics for Flury comparisons. Then load the results from 05a and compare the true statistics to the null-distribution and save the results for Table C9.

**05_OU-Process.R**
Model the trait divergence during evolution assuming a Ornstein-Uhlenbeck (OU) process. Produce Figures 8 and B9.

**07_Tensor_Analysis**
Analysis and Figures for Appendix D.

**08_Male_frequency.R**
Estimates male frequency from worm-tracking data and a trained model in 08_xgb_preds.rda for Figure B1.

**09_LD_effects_on_G.R**
**G**-matrix estimation from genomic data, varying LD weighting, for Figure 6. 

**Simulations**
This sub-directory contains all the code of the drift simulations.

**01_Adjust_G.R**
Estimate the global scaling parameter that minimizes the difference between the experimental A6140 **G**-matrix and the simulated one (see Appendix A).

**02_Main_model_code.R**
Run the simulations and process the results to obtain the final **G** matrices after 100 generations.

**03_Estimate_nb_mutation_after_evolution.R**
Simulate random neutral mutations in a population and estimate their numbers and frequencies segregating after 100 generations of drift. Produce Figure B4.

**04_Model_functions.R**
Functions used during the simulations in files 01 and 02.




