
library(zigzag)
priordat = matrix(c(-Inf, -Inf), nrow = 1)
goi = "gene1"
row.names(priordat) = c("gene1")

sq=rbind(c(1,2),c(3,4)); layout(sq)

mm <- zigzag(data = priordat, gene_length = NULL, candidate_gene_list = goi,
             output_directory = "under_prior", num_active_components =2,
             threshold_i = 0, threshold_a = c(0,5), active_gene_set = NULL,
             shared_active_variances = T, beta = 0)

mm$burnin(sample_frequency = 20, burnin_target_acceptance_rate=0.44, progress_plot = T,
          write_to_files = T, ngen=4000, append = F)


mm$mcmc(sample_frequency = 50, progress_plot = F, write_to_files = T, ngen=4000000, append = F,
        run_posterior_predictive = F, mcmcprefix = "prior1")


