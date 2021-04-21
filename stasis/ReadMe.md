#### Code and data associated with *Selection and genetic drift determine phenotypic stasis with genetic divergence*
[Bloated preprint](https://www.biorxiv.org/content/10.1101/778282v1.full).

#### 01a_Analysis_Cemee_Pop_WI.R
Load the transition rates for the CeMEE lines (`Final_merged_data_NGM.txt`), standardize data, adjusting for covariates, and estimate the **G**-matrices for each replicate population. Saves processed trait data to `Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData`.

#### 01b_Produce_Random_Gs.R
Load `Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData`, generate a null distribution of **G**-matrices by permuting line and block structure, for each replicate population (saving to `Random_**G**_Analysis_Cemee_Pop_WI_%pop_.RData`), and also for subsampling of the A6140 ancestor (saving to `Random_G_Analysis_Cemee_Pop_WI_A6140_subset.RData`).

#### 01c_Analyze_random_Gs.R
Load the observed and null **G**-matrices generated in *01b*, plot the posterior summary stats for figure 2 and figure S7.

#### 02a_State_frequency.R 
Load the trait data (`Cemee_Pop_WI/Analysis_Cemee_Pop_WI.RData`) to compare observed locomotion state frequencies to those assumed from a Markov process. Saves to `State_frequencies.RData` and makes figure S3. 

#### 02b_Compute_TR_estimates_with_freq.R 
Load `Cemee_Pop_WI/State_frequencies.RData` and the trait data (`Analysis_Cemee_Pop_WI.RData`), generates conditional mean trait values for each population sample, saves to `Cemee_Pop_WI/Analysis_Cemee_Pop_WI_betas_with_freq.RData`.

#### 02c_Compute_TR_estimates_with_freq_males.R
As above (*02c*), for males.

#### 02d_Figure1B.R
Plots figure 1B from data generated in *02b* and *02c*.

#### 02e_Figure2BC.R
Plot distributions of sequential phenotypic changes.

#### 02f_Male_frequency_and_Inbreeding.R
Fits linear models to test for effects of male frequency and inbreeding on trait values in outbred populations and inbred lines, and plots figure S8. Loads `Analysis_Cemee_Pop_WI_betas_with_freq.RData` and `population_rep_male_freqs.txt`

#### 03a_Gmatrix_Founders.R
Estimates **G** for the founders of experimental evolution. Loads `Analysis_Cemee_Pop_WI_with_gamma.RData`, saves `WI_G-matrix-all.RData`.

#### 03b_Gamma_matrix_estimation.R 
Estimate the gamma matrix from regression of fertility on conditional trait means. Plot figure 4 and figure S10. Loads `Analysis_Cemee_Pop_WI.RData`, `fertility.rda`, saves `Analysis_Cemee_Pop_WI_with_gamma.RData`.

#### 03c_Figure_5.R
Plot figure 5. Loads `Analysis_Cemee_Pop_WI_with_gamma.RData`, `WI_G-matrix-all.RData`, saves `G_matrix_projection_for_plot.RData`.

#### 04_Tensor_Analysis.R
Runs the eigentensor decomposition of observed and null **G** matrices, plots figure 3A,B and figure S11. Loads `Analysis_Cemee_Pop_WI_with_gamma.RData`, sources `Tensor_Analysis/functions_tensor.R`, saves `Cemee_Pop_WI/Analysis_Cemee_Pop_WI_tensor.RData` and intermediate data in `Tensor_Analysis/File_for_parallel_processing.RData`, `Tensor_Analysis/Tensor_processed.Rdata`.





