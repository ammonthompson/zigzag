
# set working dir to script location
library(zigzag)
library(this.path)
setwd(dirname(this.path()))

# assumes all 3 components share same variance parameter

# sim data zigzag analysis
sim_prefix <- "sim1"
out_dir    <- paste0("../data/simulated_data/zigzag_output_", sim_prefix)
simdat_fn  <- paste0("../data/simulated_data/", sim_prefix, ".tsv")
simgl_fn   <- paste0("../data/simulated_data/", sim_prefix, "_gene_length.tsv")
simdat     <- read.table(simdat_fn, header = T, row.names = 1)
simgl      <- read.table(simgl_fn, header = T, row.names = 1)

mm         <- zigzag$new(data = simdat, gene_length = simgl, output_directory = out_dir)
mm$burnin(sample_frequency = 10, ngen=5000, write_to_files = T)
mm$mcmc(sample_frequency = 10, ngen=10000, run_posterior_predictive = T, mcmcprefix = sim_prefix)


# lung data zigzag analysis
lungdat    <- read.table("../data/example_lung.tpm", header=TRUE, row.names=1)
lunggl     <- read.table("../data/example_gene_length.txt", row.names = 1, header = TRUE)

mm         <- zigzag$new(data = lungdat, gene_length = lunggl, output_directory = "../data/example_lung")
mm$burnin(sample_frequency = 10, ngen=5000, progress_plot = T)
mm$mcmc(sample_frequency = 10, ngen=10000, append = F,
        run_posterior_predictive = T, mcmcprefix = "lung")



#######################################
## mcmc under prior. (set beta = 0) ###
#######################################

priordat = matrix(c(-Inf, -Inf), nrow = 1)
row.names(priordat) = c("gene1")

mm <- zigzag$new(data = priordat,
                 output_directory = "testing_prior",
                 num_active_components = 2,
                 threshold_i = 0,
                 threshold_a = c(0,5),
                 beta = 0)

mm$burnin(sample_frequency = 20, ngen=4000, write_to_files = T)

mm$mcmc(sample_frequency = 50, ngen=4000000, mcmcprefix = "prior1")


