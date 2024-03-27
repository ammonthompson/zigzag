
# set working dir to script location
library(zigzag)
library(this.path)
setwd(dirname(this.path()))



# assumes 2 active and 1 inactive component and that all 3 components share same variance parameter
sq=rbind(c(1,2),c(3,4)); layout(sq)

# sim data zigzag analysis
simdat <- read.table("../simulate_data/simulated_data/sim1.tsv", header = T, row.names = 1)
simgl <- read.table("../simulate_data/simulated_data/sim1_gene_length.tsv", header = T, row.names = 1)

mm <- zigzag$new(data = simdat, gene_length = simgl, output_directory = "sim_test", shared_variance = F)
mm$burnin(sample_frequency = 10, write_to_files = T, ngen=2000, progress_plot = T)
mm$mcmc(sample_frequency = 10, ngen=1000, append = F,
        run_posterior_predictive = T, mcmcprefix = "sim_test")


# lung data zigzag analysis
lungdat <- read.table("../quick_start_guide/example_lung.tpm", header=TRUE, row.names=1)
lunggl = read.table("../quick_start_guide/example_gene_length.txt", row.names = 1, header = TRUE)

mm <- zigzag$new(data = lungdat, gene_length = lunggl, output_directory = "lung_test")
mm$burnin(sample_frequency = 10, ngen=2000, progress_plot = T)
mm$mcmc(sample_frequency = 10, ngen=1000, append = F,
        run_posterior_predictive = T, mcmcprefix = "lung_test")

# dmel ag data zigzag analysis
agdat <- read.table("../../../RESEARCH_PROJECTS/transcriptome_turnover/data_files/tpm_files/dmel_ag_refonly_replicates_matchGLengthfile.tpm", header=TRUE, row.names=1)
aggl = read.table("../../../RESEARCH_PROJECTS/transcriptome_turnover/data_files/gene_length_files/dmel_gene_meanLength_matchTPMfile.txt", row.names = 1, header = TRUE)

mm <- zigzag$new(data = agdat, gene_length = aggl, output_directory = "pp_ag_test", variance_g_upper_bound = 25)
mm$burnin(sample_frequency = 10, ngen=2000, progress_plot = T)
mm$mcmc(sample_frequency = 10, ngen=20000, append = F,
        run_posterior_predictive = T, mcmcprefix = "pp_ag_test")



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


