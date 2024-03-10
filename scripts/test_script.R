# set working dir to script location
library(zigzag)
library(this.path)
setwd(dirname(this.path()))

# read some data

dat <- read.table("../quick_start_guide/example_lung.tpm", header=TRUE, row.names=1)
gl = read.table("../quick_start_guide/example_gene_length.txt", row.names = 1, header = TRUE)

sq=rbind(c(1,2),c(3,4)); layout(sq)

# assumes 2 active and 1 inactive component and that all 3 components share same variance parameter
mm <- zigzag$new(data = dat[1:5000,c(1,2,3,4)], gene_length = gl[1:5000,], threshold_i = -1, shared_active_inactive_variances = F,
                 output_directory = "lung_test", num_active_components = 2, threshold_a = c(1,5))

mm$burnin(sample_frequency = 10, burnin_target_acceptance_rate=0.44,
          write_to_files = T, ngen=10000, progress_plot = T)

start_time = Sys.time()
mm$mcmc(sample_frequency = 5, write_to_files = T, ngen=10000, append = F,
        run_posterior_predictive = T, mcmcprefix = "Test")
mir_runtime = Sys.time() - start_time


mm$mirror_moves = F
start_time = Sys.time()
mm$mcmc(sample_frequency = 2, write_to_files = T, ngen=10000, append = F,
        run_posterior_predictive = F, mcmcprefix = "lung_Test_nomirror")
nomir_runtime = Sys.time() - start_time
print(mir_runtime)
print(nomir_runtime)


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




