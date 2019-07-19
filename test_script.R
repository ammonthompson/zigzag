# read some data


dat <- read.table("../simulate_data/simulated_data/sim85_experiment_lib4_var2_nactive1919.tsv", header=TRUE, row.names=1)

gene_length_df = read.table("../simulate_data/simulated_data/sim85_experiment_lib4_var2_gene_length.tsv", row.names = 1, header = FALSE)


sq=rbind(c(1,2),c(3,4)); layout(sq)

mm <- zigzag$new(data = dat, gene_length = gene_length_df, candidate_gene_list = "random",
                    output_directory = "../testing", num_active_components =2,
                    threshold_a = c(0,3), active_gene_set = NULL, shared_active_variances = T, beta = 1)


mm$burnin(sample_frequency = 20, burnin_target_acceptance_rate=0.44, progress_plot = T,
          write_to_files = T, ngen=5000, append = F)


mm$mcmc(sample_frequency = 50, progress_plot = F, write_to_files = T, ngen=100000, append = F,
        run_posterior_predictive = F, mcmcprefix = "sim85_t1")



########################################
## mcmc under prior. (set beta to 0) ###
########################################

priordat = matrix(c(-Inf, -Inf), nrow = 1)
goi = "gene1"
row.names(priordat) = c("gene1")

sq=rbind(c(1,2),c(3,4)); layout(sq)

mm <- zigzag(data = priordat, gene_length = NULL, candidate_gene_list = goi,
                    output_directory = "../testing_prior", num_active_components =2,
                    threshold_i = 0, threshold_a = c(0,5), active_gene_set = NULL,
                    shared_active_variances = T, beta = 0)

mm$burnin(sample_frequency = 20, burnin_target_acceptance_rate=0.44, progress_plot = T,
          write_to_files = T, ngen=4000, append = F)


mm$mcmc(sample_frequency = 50, progress_plot = F, write_to_files = T, ngen=4000000, append = F,
        run_posterior_predictive = F, mcmcprefix = "prior1")



