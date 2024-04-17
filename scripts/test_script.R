
# set working dir to script location
library(zigzag)
library(this.path)
setwd(dirname(this.path()))

# assumes 2 active and 1 inactive component and that all 3 components share same variance parameter

# sim data zigzag analysis
sim_prefix <- paste0("sim", i)
simdat_fn  <- paste0("../simulate_data/simulated_data/", sim_prefix, ".tsv")
simgl_fn   <- paste0("../simulate_data/simulated_data/", sim_prefix, "_gene_length.tsv")
simdat     <- read.table(simdat_fn, header = T, row.names = 1)
simgl      <- read.table(simgl_fn, header = T, row.names = 1)

mm         <- zigzag$new(data = simdat, gene_length = simgl, output_directory = sim_prefix)
mm$burnin(sample_frequency = 10, ngen=5000, write_to_files = T)
mm$mcmc(sample_frequency = 10, ngen=10000, run_posterior_predictive = T, mcmcprefix = sim_prefix)


# lung data zigzag analysis
lungdat    <- read.table("../quick_start_guide/example_lung.tpm", header=TRUE, row.names=1)
lunggl     <- read.table("../quick_start_guide/example_gene_length.txt", row.names = 1, header = TRUE)

mm         <- zigzag$new(data = lungdat, gene_length = lunggl, output_directory = "lung_test")
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
                    shared_variance_prior_min = 1, shared_variance_prior_max = 3, beta = 0)

mm$burnin(sample_frequency = 20, burnin_target_acceptance_rate=0.44, progress_plot = T,
          write_to_files = T, ngen=4000, append = F)


mm$mcmc(sample_frequency = 50,  write_to_files = T, ngen=4000000, append = F,
        run_posterior_predictive = F, mcmcprefix = "prior1")


