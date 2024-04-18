#!/usr/bin/env Rscript
library(zigzag)
args = commandArgs(trailingOnly = T)
start_idx = as.numeric(args[1])
end_idx = as.numeric(args[2])
for(i in start_idx:end_idx){
  # sim data zigzag analysis
  sim_prefix = paste0("sim", i)
  simdat_fn = paste0("../data/simulated_data/", sim_prefix, ".tsv")
  simgl_fn = paste0("../data/simulated_data/", sim_prefix, "_gene_length.tsv")
  simdat <- read.table(simdat_fn, header = T, row.names = 1)
  simgl <- read.table(simgl_fn, header = T, row.names = 1)

  mm <- zigzag$new(data = simdat, gene_length = simgl, output_directory = sim_prefix, shared_variance = T,
                   threshold_a = c(1, 6), threshold_i = -1, num_active_components = 2,
                   weight_active_shape_1 = 10, weight_active_shape_2 = 10,
                   spike_prior_shape_1 = 2, spike_prior_shape_2 = 8,
                   tau_shape = 5, tau_rate = 20, s0_mu = -1, s0_sigma = 1/5, s1_shape = 1, s1_rate = 10,
                   alpha_r_shape = 15, alpha_r_rate = 1, inactive_means_prior_shape =1, inactive_means_prior_rate = 1,
                   active_means_dif_prior_shape = 1, active_means_dif_prior_rate = 1,
                   shared_variance_prior_min = 1, shared_variance_prior_max = 3)
  mm$burnin(sample_frequency = 10, write_to_files = T, ngen=5000, progress_plot = F)
  mm$mcmc(sample_frequency = 10, ngen=20000, append = F,
          run_posterior_predictive = F, mcmcprefix = sim_prefix)
}

