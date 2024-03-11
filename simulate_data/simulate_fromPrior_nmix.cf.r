#INPUTS:
#	1:	control file
#	2:	Output file prefix
#	3:	How many randomly selected Yg values to output
##########################

file_prefix = "sim"
num_yg = 500
df = read.table("sim_control_file.template.cf",header = F, row.names = 1)

simulated_gene_lengths  = rlnorm(df$V2[5], log(2000), 0.5)

output_dir = "simulated_data/"
params = as.numeric(df[,1])

num_libs = params[4]

#### Functions
get_sigma2_g = function(gene){

    sigma_x = exp(s0x + s1x * gene)

    return(rlnorm(length(gene), log(sigma_x) + taux, sqrt(taux)))

}


get_p_x = function(gene, lib, gl){

  return(1 - exp(-alpha_rx[lib] * gl * exp(gene)))

}

prior_settings_file = paste0(output_dir, file_prefix, ".parameterValues")
write.table(matrix(c("sim_number", "mu_i", "variance_i", "spike_prob", "weight_active",
	"mu_a1", "mu_a2", "variance_a1", "variance_a2", "w_a1", "w_a2", "tau", "s0", "s1",
	paste0(rep("alpha_r", num_libs), seq(num_libs))), nrow = 1),
	file = prior_settings_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

num_transcripts = params[5]
gene_names = paste0("gene_", seq(num_transcripts))

### Generate sims
for(i in seq(params[17])){



	#### spike and weights
	total = params[5]
	sim_gl = sample(simulated_gene_lengths, total)
	sim_rel_gl = sim_gl/mean(sim_gl)

	weight_active =  params[13]
	if(is.na(params[13])) weight_active =  rbeta(1, 10, 10)

	weight_within_active = params[14:15]
	if(is.na(params[14]) | is.na(params[15])){weight_within_active[1] = rbeta(1, 10, 2); weight_within_active[2] = 1 - weight_within_active[1]}

	spike_prob =  params[16]
	if(is.na(params[16])) spike_prob = rbeta(1, 2, 8)

	num_inactive = floor((1 - weight_active) * total) + rbinom(1, 1, (1 - weight_active) * total - floor((1 - weight_active) * total))

        num_spike = floor(spike_prob * num_inactive) + rbinom(1, 1, spike_prob * num_inactive - floor(spike_prob * num_inactive))

	#### sigma_x
	taux = params[1]; if(is.na(taux)) taux = rgamma(1, 4, 40)
	s0x = params[2]; if(is.na(s0x)) s0x = rnorm(1, -3, 1/5)
	s1x = params[3]; if(is.na(s1x)) s1x = -rgamma(1, 1, 10)


	#### p_x

	alpha_rx = rep(params[6], num_libs); if(is.na(params[6])) alpha_rx = rgamma(num_libs, 15, 1)


	### sim_yg

	ymu_i = params[7]; if(is.na(ymu_i)) ymu_i = -rgamma(1, 1, 1)
	ymu_a1 = params[8]; if(is.na(ymu_a1)) ymu_a1 = 1 + rgamma(1, 1, 1)
	ymu_a2 = params[9]; if(is.na(ymu_a2)) ymu_a2 = 7 + rgamma(1, 1, 1)

	yVariance_i =  params[10]; if(is.na(yVariance_i)) yVariance_i = exp(runif(1, log(1), log(3)))
	yVariance_a1 =  params[11]; if(is.na(yVariance_a1)) yVariance_a1 = exp(runif(1, log(1), log(3)))
	yVariance_a2 =  params[12]; if(is.na(yVariance_a2)) yVariance_a2 = yVariance_a1

	cat("simulation:", i, "\n")
	cat("ymu_i:", round(ymu_i, digits = 5), " yVariance_i:", round(yVariance_i, digits = 5), " spike prob:", round(spike_prob, digits = 5), " weight_active:", round(weight_active, digits = 5))
	cat(" ymu_a1:", round(ymu_a1, digits = 5), " yVariance_a1:", round(yVariance_a1, digits = 5), " weight_within_active:", round(weight_within_active, digits = 5))
	cat(" ymu_a2:", round(ymu_a2, digits = 5), " yVariance_a2:", round(yVariance_a2, digits = 5))
	cat(" taux:", round(taux, digits = 5), " s0x:", round(s0x, digits = 5), " s1x:", round(s1x, digits = 5), " alpha_r:", round(alpha_rx, digits = 5),  "\n")
	print("  ")


	write.table(matrix(c(paste0("sim_", i), round(c(ymu_i, yVariance_i, spike_prob, weight_active, ymu_a1, ymu_a2, yVariance_a1, yVariance_a2, weight_within_active, taux, s0x, s1x, alpha_rx), digits = 3)), nrow = 1),
	file = prior_settings_file, append = TRUE, sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

	num_a1 = round((total - num_inactive) * weight_within_active[1], 0)
	num_a2 = (total - num_inactive) - num_a1


	### simulate Yg
	sim_yg = c(rnorm(num_inactive - num_spike, ymu_i, sqrt(yVariance_i)), rnorm(num_a1, ymu_a1, sqrt(yVariance_a1)), rnorm(num_a2, ymu_a2, sqrt(yVariance_a2)))

	### simulate Xg
	sigma_gx <- get_sigma2_g(sim_yg)

	sim_xg = sapply(1:num_libs, function(x){return(rnorm(length(sim_yg), sim_yg, sqrt(sigma_gx)))})
	#one lib is an outlier for many genes
	# sim_xg[sample(seq(nrow(sim_xg)), 2000, replace = F),1] = rnorm(2000, 4, 2)

	sapply(1:num_libs, function(lib){

	  p_x = get_p_x(sim_yg, lib, sim_rel_gl[(num_spike + 1):total])

	  sim_xg[(runif(length(sim_yg)) < (1 - p_x)), lib] <<- -Inf  #because the in-spike genes will be place at the beginning of sim_xg/yg

	})


	### simulate spike # first num_spike genes are in spike
	sim_xg = as.data.frame(rbind( matrix(rep(rep(-Inf, num_libs),num_spike), nrow=num_spike), sim_xg))

        row.names(sim_xg) = gene_names

	sim_yg=c( rep(-Inf,num_spike), sim_yg)
	sigma_gx = c(rep(1, num_spike), sigma_gx)

	nactive = num_a1 + num_a2

	colnames = c("gene_id", sapply(seq(num_libs), function(x){return(paste("lib", x, sep ="_"))}))


	# sample random genes to record Yg values. Store in list for zigzag input as goi
	goi = gene_names[sample(seq(num_transcripts), num_yg, replace = FALSE)]
	write.table(gene_names[which(gene_names %in% goi)], file = paste0(output_dir, "sim", i,"_goi.txt"), row.names = FALSE,
	            col.names = FALSE, sep = "\t", quote = FALSE)


	# expression data file
	write.table(cbind(row.names(sim_xg), exp(sim_xg)),
	            file = paste0(output_dir, "sim", i, "_", file_prefix, "_nactive", nactive, ".tsv"),
	            quote = F, row.names = FALSE, col.names = colnames, sep="\t")

	# gene length file
  write.table(cbind(row.names(sim_xg), sim_gl),
              file = paste0(output_dir, "sim", i, "_", file_prefix, "_gene_length.tsv"),
              quote=FALSE, row.names = FALSE, col.names=c("gene_name", "gene_length"), sep = "\t")

  # yg and simga_g of interest data
	write.table(cbind("gene", "Yg", "sigma_g"),
	            file = paste0(output_dir, "sim", i,"_Yg.txt"),
	            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

	write.table(cbind(gene_names[which(gene_names %in% goi)],
	                  sim_yg[which(gene_names %in% goi)],
	                  sigma_gx[which(gene_names %in% goi)]),
	            file = paste0(output_dir, "sim", i,"_Yg.txt"),
	            row.names = FALSE,	quote = FALSE, append = TRUE, sep = "\t", col.names = FALSE )


#	print(warnings())
#	warnings()
}
warnings()
print(warnings())
