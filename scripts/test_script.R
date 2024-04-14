
# set working dir to script location
library(zigzag)
library(this.path)
setwd(dirname(this.path()))



# assumes 2 active and 1 inactive component and that all 3 components share same variance parameter
sq=rbind(c(1,2),c(3,4)); layout(sq)


for(i in 37:69){
  # sim data zigzag analysis
  sim_prefix = paste0("sim", i)
  simdat_fn = paste0("../simulate_data/simulated_data/", sim_prefix, ".tsv")
  simgl_fn = paste0("../simulate_data/simulated_data/", sim_prefix, "_gene_length.tsv")
  simdat <- read.table(simdat_fn, header = T, row.names = 1)
  simgl <- read.table(simgl_fn, header = T, row.names = 1)

  mm <- zigzag$new(data = simdat, gene_length = simgl, output_directory = "check_sim", shared_variance = T,
                   threshold_a = c(1, 7), threshold_i = -1, num_active_components = 2,
                   weight_active_shape_1 = 10, weight_active_shape_2 = 10,
                   spike_prior_shape_1 = 2, spike_prior_shape_2 = 8,
                   tau_shape = 4, tau_rate = 20, s0_mu = -2, s0_sigma = 1/5, s1_shape = 1, s1_rate = 10,
                   alpha_r_shape = 15, alpha_r_rate = 1, inactive_means_prior_shape =1, inactive_means_prior_rate = 1,
                   active_means_dif_prior_shape = 1, active_means_dif_prior_rate = 1,
                   shared_variance_prior_min = 1, shared_variance_prior_max = 3)
  mm$burnin(sample_frequency = 10, write_to_files = T, ngen=5000, progress_plot = F)
  mm$mcmc(sample_frequency = 10, ngen=20000, append = F,
          run_posterior_predictive = F, mcmcprefix = sim_prefix)
}

# lung data zigzag analysis
lungdat <- read.table("../quick_start_guide/example_lung.tpm", header=TRUE, row.names=1)
lunggl = read.table("../quick_start_guide/example_gene_length.txt", row.names = 1, header = TRUE)

mm <- zigzag$new(data = lungdat, gene_length = lunggl, output_directory = "lung_test")
mm$burnin(sample_frequency = 10, ngen=2000, progress_plot = T)
mm$mcmc(sample_frequency = 10, ngen=1000, append = F,
        run_posterior_predictive = T, mcmcprefix = "lung_test")



########################################
## mcmc under prior. (set beta to 0) ###
########################################

priordat = matrix(c(-Inf, -Inf), nrow = 1)
goi = "gene1"
row.names(priordat) = c("gene1")

sq=rbind(c(1,2),c(3,4)); layout(sq)

mm <- zigzag$new(data = priordat, gene_length = NULL, candidate_gene_list = goi,
                    output_directory = "testing_prior", num_active_components =2,
                    threshold_i = 0, threshold_a = c(0,5), active_gene_set = NULL,
                    shared_active_variances = T, beta = 0, library_bias = T, bias_var = 0.001)

mm$burnin(sample_frequency = 20, burnin_target_acceptance_rate=0.44, progress_plot = T,
          write_to_files = T, ngen=4000, append = F)


mm$mcmc(sample_frequency = 50,  write_to_files = T, ngen=4000000, append = F,
        run_posterior_predictive = F, mcmcprefix = "prior1")


