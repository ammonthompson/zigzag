# read some data


# dat <- read.table("../../../bitbucket_repos/nmix_data_files/gtex_data/allgenes_liver_gtex.tpm", header=TRUE, row.names=1)
dat <- read.table("../simulate_data/simulated_data/sim3_twoComps_2mu3_gene_length.tsv", header=TRUE, row.names=1)

gene_length_df = read.table("../simulate_data/simulated_data/sim3_twoComps_2mu3_gene_length.tsv", row.names = 1, header = TRUE)


sq=rbind(c(1,2),c(3,4)); layout(sq)

mm <- zigzag$new(data = dat, gene_length = gene_length_df, candidate_gene_list = "random",
                    output_directory = "../testing", num_active_components =2,
                    active_means_dif_prior_shape = 20000, active_means_dif_prior_rate = 10000,
                    threshold_i = -1, threshold_a = c(0, 1), active_gene_set = NULL, shared_active_variances = T, beta = 1)


mm$burnin(sample_frequency = 20, burnin_target_acceptance_rate=0.44, progress_plot = T,
          write_to_files = T, ngen=5000, append = F)


mm$mcmc(sample_frequency = 50, progress_plot = F, write_to_files = T, ngen=50000, append = F,
        run_posterior_predictive = T, mcmcprefix = "sim31")



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



