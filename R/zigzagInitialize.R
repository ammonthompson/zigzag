zigzag$methods(

  ###############
  # Constructor #
  ###############

  initialize = function(data, gene_length = NULL, candidate_gene_list = "all",
                        num_active_components = "auto",
                        weight_active_shape_1 = 2,
                        weight_active_shape_2 = 2,
                        shared_variance_prior_min = 0.01,
                        shared_variance_prior_max = 5,
                        inactive_means_prior_shape = 1,
                        inactive_means_prior_rate = 1/3,
                        inactive_variances_prior_min = shared_variance_prior_min,
                        inactive_variances_prior_max = shared_variance_prior_max,
                        spike_prior_shape_1 = 1,
                        spike_prior_shape_2 = 1,
                        active_means_dif_prior_shape = 1,
                        active_means_dif_prior_rate = 1/3,
                        active_variances_prior_min = shared_variance_prior_min,
                        active_variances_prior_max = shared_variance_prior_max,
                        shared_active_variance = TRUE,
                        shared_variance = TRUE,
                        output_directory = "output",
                        multi_ta = FALSE,
                        threshold_a = "auto",
                        threshold_i = "auto",
                        beta = 1,
                        temperature = 1,
                        tau_shape = 1,
                        tau_rate = 1,
                        s0_mu = -1,
                        s0_sigma = 2,
                        s1_shape = 1,
                        s1_rate = 2,
                        variance_g_upper_bound = Inf,
                        alpha_r_shape = 1,
                        alpha_r_rate = 1/10,
                        active_gene_set = NULL,
                        ...) {

    .self$shared_variance <- shared_variance
    .self$shared_active_variance <- shared_active_variance
    if(shared_variance) .self$shared_active_variance <- TRUE


    ############################
    ## Set up data matrix, Xg ##
    ############################

    if(min(data) >= 0) .self$Xg <- log(as.matrix(data)) else .self$Xg <- as.matrix(data)

    .self$inf_tol <- -Inf

    .self$Xg[Xg == -Inf] <- inf_tol


    ########################################################
    ## Set thresholds for component mean priors ############
    ########################################################

    # if any of the component mean priors are set to "auto" then auto find prior thresholds
    # use for each that apply
    if(any(is.character(num_active_components), is.character(threshold_i), is.character(threshold_a))){
      rmedians_noInf <- matrixStats::rowMedians(Xg)
      rmedians_noInf <- rmedians_noInf[rmedians_noInf > -Inf]
      auto_thresholds_list <- .self$findthresholds(rmedians_noInf, num_active_components)
    }

    # set inactive mean prior threshold
    if(any(is.character(threshold_i), threshold_i == "auto")){
      thresh_i <- auto_thresholds_list[[2]][1]
    }else if(length(threshold_i) == 1){
      thresh_i <- threshold_i
    }else{
      cat("\tWarning: Currently, only 1 inactive component allowed. Threshold set to first in vector.\n")
      thresh_i <- threshold_i[1]
    }
    .self$threshold_i <- thresh_i

    # set number of active components.
    # if "auto" or character, get from auto_thresholds_list
    if(is.character(num_active_components)){
      num_acomps <- auto_thresholds_list[[1]]
    }else{
      num_acomps <- num_active_components
    }
    .self$num_active_components <- num_acomps

    # set active mean prior threshold
    # if auto or character, get from auto_thresholds_list
    if(any(is.character(threshold_a), threshold_a == "auto")){
      thresh_a <- auto_thresholds_list[[2]]
    }else{
      thresh_a <- threshold_a
    }
    .self$threshold_a <- thresh_a

    # if the length of threshold_a set by user (implied number of active components)
    # and number of active components disagree, then set the number of components equal to the length
    # of threshold_a and print a warning.
    if(length(.self$threshold_a) != .self$num_active_components){
      cat("\tWarning: num. active components does not equal length of threshold_a.
                 Number of active components set to length of threshold_a.\n")
      num_acomps <- length(threshold_a)
      .self$num_active_components <- num_acomps
    }


    if(length(thresh_a) == num_acomps){
      multi_ta <- TRUE
      .self$multi_thresh_a <- multi_ta
    }else{
      .self$multi_thresh_a <- multi_ta
    }



    #################################
    ## Set up general data params  ##
    #################################

    cat("Loading and log-transforming data...\n")

    .self$sqrt2pi <- sqrt(2*pi)
    .self$temperature <- temperature
    .self$beta <- beta
    .self$beta_oneLib <- 1
    .self$gen <- 0
    .self$gene_names <- rownames(data)
    .self$num_libraries <- ncol(data)
    .self$num_transcripts <- nrow(data)

    if(num_libraries < 2) stop("Data must have > 1 library.")

    # components matrix
    .self$component_matrix <- matrix(rep(c(0,seq(num_acomps)), num_transcripts),
                                     ncol=num_acomps+1, byrow = T)



    # Gene lengths, scaled to have mean = 1
    if(! is.null(gene_length)){
      if( is.vector(gene_length) ) .self$gene_lengths <- gene_length/mean(gene_length) else .self$gene_lengths <- gene_length[,1]/mean(gene_length[,1])
    }else{
      .self$gene_lengths <- rep(1, num_transcripts)
    }
    if(length(gene_lengths) != num_transcripts) cat("WARNING: genelengths size does not equal number of genes!!!\n")


    # candidate genes to output posterior samples to files
    cat("Reading candidate genes to record in mcmc log files...\n")
    if(candidate_gene_list[1] == "random"){
      .self$candidate_gene_list <- c(gene_names[c(1:10,sample(num_transcripts, size = 100, replace = F))])
    }else if(candidate_gene_list[1] == "all"){
      .self$candidate_gene_list <- gene_names
    }else{
      .self$candidate_gene_list <- candidate_gene_list
      if(length(which(!( candidate_gene_list %in% gene_names))) > 0){
        stop("Invalid candidate_gene_list. Check for typos and gene names not present in the data.")
      }
    }

    # Set up active genes: the genes user is certain are actively expressed. Default is NULL
    if(!is.null(active_gene_set)){

      .self$active_gene_set <- active_gene_set
      .self$active_gene_set_idx <- sapply(seq(active_gene_set), function(gene) which(gene_names == active_gene_set[gene]) )

    }

    # create output directory if doesnt exist
    .self$output_directory <- output_directory
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = TRUE, showWarnings = TRUE)
      cat("Created output directory called: ", output_directory, "...\n")
    }


    ##########################################################
    ### initialize proposal tuning parameters "tuningParam" ##
    ##########################################################

    .self$inactive_mean_tuningParam <- 0.5
    .self$inactive_variance_tuningParam <- 0.5
    .self$spike_probability_tuningParam <- 1000
    .self$active_mean_tuningParam <- rep(0.5, num_acomps)
    .self$active_variance_tuningParam <- rep(0.5, num_acomps)
    .self$mixture_weight_tuningParam <- 1000
    .self$tuningParam_s0 <- 0.05
    .self$tuningParam_s1 <- 0.05
    .self$tuningParam_tau <- 0.05
    .self$tuningParam_s0tau <- 0.05
    .self$tuningParam_alpha_r <- rep(0.5, num_libraries)
    .self$tuningParam_yg <- rep(0.5, num_transcripts)
    .self$tuningParam_variance_g <- rep(0.5, num_transcripts)
    .self$tuningParam_multi_sigma <- 0.5
    .self$tuningParam_sigma_mu <- 0.5



    #############################################################
    ## Set up parameter proposal relative probabilities #########
    #############################################################

    .self$proposal_list <- list(.self$gibbsMixtureWeights,
                                .self$gibbsAllocationActiveInactive,
                                .self$gibbsAllocationWithinActive,
                                .self$mhInactiveMeans,
                                .self$mhInactiveVariances,
                                .self$mhActiveMeansDif,
                                .self$mhActiveVariances,
                                .self$gibbsSpikeProb,
                                .self$mhSpikeAllocation,
                                .self$mhYg,
                                .self$mhVariance_g,
                                .self$mhTau,
                                .self$mhSg,
                                .self$mhS0Tau,
                                .self$mhP_x)

    if(shared_variance) .self$proposal_list[[7]] <- .self$mhSharedVariance

    is2Libs <- (num_libraries == 2) * 0.5 # use this variable to upweight probability of proposing L1 variance params

    if(num_libraries > 1 ){

      .self$proposal_probs <- c(4, 30, 5,                                      ### weights, alloc active_inactive, alloc within_active
                                6, 7 * (1 - shared_variance),                      ### i_mean, i_var
                                5,                                                  ### a_mean
                                5 + 2 * (1 - shared_active_variance),                ### a_var
                                2, 5,                                               ### spike prob, spike alloc
                                c(1, 1) + is2Libs + 1 * (num_transcripts < 15000),   ### Yg, sigm_g
                                c(8,10,8),                                        ### tau, Sg, s0tau
                                num_libraries * 1.1                                  ### p_x
      )

    }else{
      # not implemented yet
      .self$proposal_probs <- c(8, 60,10,                                                      ### weights, alloc active_inactive, alloc within_active
                                12, 15 * (1 - shared_variance),               ### i_mean, i_var
                                8,                                                             ### a_mean
                                4 + 4 * (1 - shared_active_variance),                         ### a_var
                                5, 0,                                                          ### spike prob, spike alloc
                                c(1, 1) * 0,                                                   ### Yg, sigm_g
                                c(6, 6, 6) * 0,                                                ### tau, Sg, s0tau
                                num_libraries * 0.75 * 0)                                      ### p_x

    }


    ##########################
    # Initialize the priors. #
    # Level 2 (upper) priors #
    ##########################

    cat("Initializing upper level parameters...\n")

    # active vs inactive
    .self$weight_active_shape_1 <- weight_active_shape_1
    .self$weight_active_shape_2 <- weight_active_shape_2
    .self$weight_active <- rbeta(1, weight_active_shape_1 + 1, weight_active_shape_2 + 1)
    .self$weight_active_proposed <- weight_active

    .self$weight_active_prob <- .self$computeActiveWeightPriorProbability(weight_active)
    .self$weight_active_prob_proposed <- .self$computeActiveWeightPriorProbability(weight_active_proposed)


    # spike priors
    .self$spike_prior_shape_1 <- spike_prior_shape_1
    .self$spike_prior_shape_2 <- spike_prior_shape_2
    .self$spike_probability <- rbeta(1, spike_prior_shape_1, spike_prior_shape_2)
    .self$spike_probability_proposed <- spike_probability
    .self$spike_probability_prob <- .self$computeSpikePriorProbability(spike_probability)
    .self$spike_probability_prob_proposed <- .self$computeSpikePriorProbability(spike_probability_proposed)
    .self$spike_probability_trace[[1]] <- lapply(1,function(x){return(c(rep(0,77),rep(1,23)))})

    # active means
    .self$active_means_dif_prior_shape <- active_means_dif_prior_shape
    .self$active_means_dif_prior_rate <- active_means_dif_prior_rate
    .self$active_means_dif <- rgamma(num_acomps, shape = active_means_dif_prior_rate, rate = active_means_dif_prior_rate)

    .self$active_means_dif_proposed <- active_means_dif
    .self$active_means_dif_prob <- .self$computeActiveMeansDifPriorProbability(active_means_dif)
    .self$active_means_dif_prob_proposed <- .self$computeActiveMeansDifPriorProbability(active_means_dif_proposed)
    .self$active_means_trace[[1]] <- lapply(1,function(x){return(t(sapply(1:num_acomps,function(y){return(c(rep(0,77),rep(1,23)))})))})
    .self$active_means <- .self$calculate_active_means(active_means_dif, mt = multi_ta)

    # shared variances
    .self$shared_variance_prior_min <- shared_variance_prior_min
    .self$shared_variance_prior_max <- shared_variance_prior_max
    .self$shared_variance_prior_log_min <- log(shared_variance_prior_min, 10)
    .self$shared_variance_prior_log_max <- log(shared_variance_prior_max, 10)

    # active variances
    if(shared_variance){
      .self$active_variances_prior_min <- shared_variance_prior_min
      .self$active_variances_prior_max <- shared_variance_prior_max
      .self$active_variances_prior_log_min <- shared_variance_prior_log_min
      .self$active_variances_prior_log_max <- shared_variance_prior_log_max
      .self$active_variances <- rep(10^(runif(1, active_variances_prior_log_min, active_variances_prior_log_max)), num_acomps)
      .self$active_variances_proposed <- active_variances
      .self$active_variances_prob <- .self$computeActiveVariancesPriorProbability(active_variances[1])
      .self$active_variances_prob_proposed <- .self$computeActiveVariancesPriorProbability(active_variances_proposed)
      .self$active_variances_trace[[1]] <- lapply(1,function(x){return(t(sapply(1:num_acomps,function(y){return(c(rep(0,77),rep(1,23)))})))})
    }else{
      .self$active_variances_prior_min <- active_variances_prior_min
      .self$active_variances_prior_max <- active_variances_prior_max
      .self$active_variances_prior_log_min <- log(active_variances_prior_min, 10)
      .self$active_variances_prior_log_max <- log(active_variances_prior_max, 10)
      if(.self$shared_active_variance){
        .self$active_variances <- rep(10^(runif(1, active_variances_prior_log_min, active_variances_prior_log_max)), num_acomps)
        .self$active_variances_proposed <- active_variances
        .self$active_variances_prob <- .self$computeActiveVariancesPriorProbability(active_variances[1])
        .self$active_variances_prob_proposed <- .self$computeActiveVariancesPriorProbability(active_variances_proposed)
        .self$active_variances_trace[[1]] <- lapply(1,function(x){return(t(sapply(1:num_acomps,function(y){return(c(rep(0,77),rep(1,23)))})))})
      }else{
        .self$active_variances <- 10^(runif(num_acomps, active_variances_prior_log_min, active_variances_prior_log_max))
        .self$active_variances_proposed <- active_variances
        .self$active_variances_prob <- .self$computeActiveVariancesPriorProbability(active_variances)
        .self$active_variances_prob_proposed <- .self$computeActiveVariancesPriorProbability(active_variances_proposed)
        .self$active_variances_trace[[1]] <- lapply(1,function(x){return(t(sapply(1:num_acomps,function(y){return(c(rep(0,77),rep(1,23)))})))})
      }
    }

    # inactive means
    .self$inactive_means_prior_shape <- inactive_means_prior_shape
    .self$inactive_means_prior_rate <- inactive_means_prior_rate
    .self$inactive_means <- thresh_i - rgamma(1, shape = 1, rate = 1)

    .self$inactive_means_proposed <- inactive_means
    .self$inactive_means_prob <- .self$computeInactiveMeansPriorProbability(inactive_means)
    .self$inactive_means_prob_proposed <- .self$computeInactiveMeansPriorProbability(inactive_means_proposed)
    .self$inactive_means_trace[[1]] <- lapply(1,function(x){return(c(rep(0,77),rep(1,23)))})

    # inactive variances
    if(shared_variance){
      .self$inactive_variances_prior_min <- shared_variance_prior_min
      .self$inactive_variances_prior_max <- shared_variance_prior_max
      .self$inactive_variances_prior_log_min <- shared_variance_prior_log_min
      .self$inactive_variances_prior_log_max <- shared_variance_prior_log_max
      .self$inactive_variances <- active_variances[1]
    }else{
      .self$inactive_variances_prior_min <- inactive_variances_prior_min
      .self$inactive_variances_prior_max <- inactive_variances_prior_max
      .self$inactive_variances_prior_log_min <- log(inactive_variances_prior_min, 10)
      .self$inactive_variances_prior_log_max <- log(inactive_variances_prior_max, 10)
      .self$inactive_variances <- 10^(runif(1, inactive_variances_prior_log_min, inactive_variances_prior_log_max))
    }
    .self$inactive_variances_proposed <- inactive_variances
    .self$inactive_variances_prob <- .self$computeInactiveVariancesPriorProbability(inactive_variances)
    .self$inactive_variances_prob_proposed <- .self$computeInactiveVariancesPriorProbability(inactive_variances_proposed)
    .self$inactive_variances_trace[[1]] <- lapply(1,function(x){return(c(rep(0,77),rep(1,23)))})


    # mixture weights within active component
    .self$weight_within_active_alpha <- rep(weight_active_shape_1/num_acomps, num_acomps)

    .self$weight_within_active <- .self$r_dirichlet(.self$weight_within_active_alpha)
    .self$weight_within_active_proposed <- weight_within_active
    .self$weight_within_active_prob <- .self$computeWithinActiveWeightPriorProbability(weight_within_active)
    .self$weight_within_active_prob_proposed <- .self$computeWithinActiveWeightPriorProbability(weight_within_active_proposed)
    .self$mixture_weight_trace[[1]] <- list(c(rep(0,77),rep(1,23)))

    # initialize the allocations
    .self$allocation_active_inactive <- rbinom(num_transcripts, size=1, p=weight_active)

    if(!is.null(active_gene_set)) .self$allocation_active_inactive[active_gene_set_idx] <- as.integer(1)
    .self$allocation_active_inactive_proposed <- allocation_active_inactive

    .self$allocation_within_active[[1]] <- sample.int(num_acomps, size=num_transcripts, replace=TRUE, prob=weight_within_active)
    .self$allocation_within_active_proposed <- allocation_within_active

    all_allocation = allocation_active_inactive * allocation_within_active[[1]]

    .self$allocation_trace <- matrix(apply(component_matrix, 2, function(comp_matrix_col) 1 *
                                             (comp_matrix_col == all_allocation)), nrow = num_transcripts)

    # combined move parameters
    .self$multi_sigma_trace[[1]] <- list(c(rep(0,77),rep(1,23)))
    .self$sigma_mu_trace[[1]] <- list(c(rep(0,77),rep(1,23)))


    ##########################
    # level 1 (lower) priors #
    ##########################

    cat("Initializing lower level parameters...\n")

    # gene detection prob
    .self$alpha_r_shape <- alpha_r_shape
    .self$alpha_r_rate <- alpha_r_rate
    .self$alpha_r <- rgamma(num_libraries, 1, 1)
    .self$alpha_r_trace[[1]] <- lapply(1,function(x){return(t(sapply(1:num_libraries,function(y){return(c(rep(0,77),rep(1,23)))})))})
    alpha_r_max = max(alpha_r)

    # gene-wise variance
    .self$s0_mu <- s0_mu
    .self$s0_sigma <- s0_sigma
    .self$s1_rate <- s1_rate
    .self$s1_shape <- s1_shape
    .self$tau_shape <- tau_shape
    .self$tau_rate <- tau_rate

    .self$s0 <- rnorm(1, 0, 0.1)
    .self$s1 <- -rgamma(1, 1, 5)
    .self$tau <- rgamma(1, 2, 5)

    .self$s0_trace[[1]] <- lapply(1,function(x){return(c(rep(0,77),rep(1,23)))})
    .self$s1_trace[[1]] <- lapply(1,function(x){return(c(rep(0,77),rep(1,23)))})
    .self$tau_trace[[1]] <- lapply(1,function(x){return(c(rep(0,77),rep(1,23)))})
    .self$s0tau_trace[[1]] <- lapply(1, function(x){return(c(rep(0,77), rep(1,23)))})

    ## Yg initialize
    .self$no_detect_spike <- ((rowSums(as.matrix((Xg[seq(num_transcripts),, drop = F] ==
                                                    rep(-Inf, num_libraries)) * 1)) == num_libraries) * 1)
    .self$inactive_spike_allocation <- no_detect_spike
    .self$in_spike_idx <- which(no_detect_spike == 1)
    .self$out_spike_idx <- which(no_detect_spike == 0)
    .self$all_zero_idx <- which(no_detect_spike == 1)
    .self$allocation_active_inactive[in_spike_idx] <- as.integer(0)
    .self$active_idx <- which(allocation_active_inactive == 1)
    .self$inactive_idx <- which(allocation_active_inactive == 0)


    # initialize Yg and variance_g near observed values in data. If undetected, sample
    if(num_libraries > 1 & max(Xg) > -Inf){

      #Yg from Xg means
      .self$rwm <- sapply(seq(num_transcripts), function(x){
        if(matrixStats::count(Xg[x,], value = -Inf) < num_libraries){
          v = Xg[x,]
          mean(v[v > -Inf])
        }else{
          rnorm(1, inactive_means, sqrt(inactive_variances))
        }
      })

      #Yg from Xg variances
      .self$rwv <- sapply(seq(num_transcripts), function(x){
        if(matrixStats::count(Xg[x,], value = -Inf) < (num_libraries - 1)){
          v = Xg[x,]
          var(v[v > -Inf])
        }else{
          rlnorm(1, log(mean(matrixStats::rowVars(Xg), na.rm = T)), sqrt(var(matrixStats::rowVars(Xg), na.rm = T)))
        }
      })
      .self$Yg <- rnorm(num_transcripts, rwm, sqrt(rwv))

    }else{

      .self$Yg <- Xg[,1]
      .self$Yg[which(Yg == inf_tol)] <- rnorm(length(which(Yg == inf_tol)), inactive_means, sqrt(inactive_variances))

    }

    .self$Yg_proposed <- Yg

    .self$Yg_trace <- t(sapply(seq(num_transcripts), function(g){return(c(rep(0,77),rep(1,23)))}))


    ## Initialize Sigma_x and Variance_g gene variance and shrinkage prior parameters
    .self$Sg <- exp(s0 + s1 * Yg)
    .self$variance_g <- rlnorm(num_transcripts, 1/2, 1/5)
    .self$variance_g_trace <- t(sapply(seq(num_transcripts), function(g){return(c(rep(0,77),rep(1,23)))}))
    .self$variance_g_upper_bound <- variance_g_upper_bound
    .self$variance_g[which(variance_g > variance_g_upper_bound)] <- variance_g_upper_bound

    .self$p_x <- .self$get_px()

    # Initialize XgLikelihood. resample Yg and downstream values if -Inf likelihood results
    .self$XgLikelihood <- .self$computeXgLikelihood(Xg, Yg, variance_g, p_x)

    cat("Reinitialize genes with zero likelihood: ", which(XgLikelihood == -Inf), "\n")

    # Occasionally there are combinations of intial values of Yg, variance_g and alpha_r
    # that can lead to 0 likelihood for some genes.
    # For those genes, reinitialize until a combination that has > 0 likelihood is generated.
    if(num_libraries > 1){

      counter = 0

      while(matrixStats::count(XgLikelihood, value = -Inf) > 0){

        .self$Yg[which(XgLikelihood == -Inf)] <- rnorm(length(which(XgLikelihood == -Inf)),
                                                       inactive_means, sqrt(inactive_variances))
        .self$Sg <- exp(s0 + s1 * Yg)
        .self$p_x <- .self$get_px()
        .self$variance_g[which(XgLikelihood == -Inf)] <- rlnorm(length(which(XgLikelihood == -Inf)),
                                                                log(Sg[which(XgLikelihood == -Inf)]) + tau, sqrt(tau))
        .self$variance_g[which(variance_g > variance_g_upper_bound)] <- variance_g_upper_bound
        .self$XgLikelihood <- .self$computeXgLikelihood(Xg, Yg, variance_g, p_x)

        counter  = counter + 1

        if(counter > 1000){
          cat("initialization failed. Try again\n")
          break
        }
      }
    }

    .self$variance_g_probability <- .self$computeVarianceGPriorProbability(variance_g, Sg)



    ## write user specified prior parameter values to file
    if(shared_variance){
      out_prior_settigs <- data.frame(cbind(c("s0_mu", "s0_sigma", "s1_shape", "s1_rate", "tau_rate", "tau_shape", "alpha_r_shape", "alpha_r_rate",
                                              "weight_active_shape_1", "weight_active_shape_2", "weight_within_active_alpha", "spike_prior_shape_1", "spike_prior_shape_2",
                                              "active_means_dif_prior_shape", "active_meaans_dif_prior_rate",
                                              "inactive_means_prior_shape", "inactive_means_prior_rate", "shared_variance_prior_min", "shared_variance_prior_max",
                                              "threshold_i", "threshold_a"),
                                            c(.self$s0_mu, .self$s0_sigma, .self$s1_shape, .self$s1_rate, .self$tau_rate, .self$tau_shape,  .self$alpha_r_shape, .self$alpha_r_rate,
                                              .self$weight_active_shape_1, .self$weight_active_shape_2, .self$weight_within_active_alpha[1], .self$spike_prior_shape_1, .self$spike_prior_shape_2,
                                              .self$active_means_dif_prior_shape, round(.self$active_means_dif_prior_rate, digits = 3) ,
                                              .self$inactive_means_prior_shape, round(.self$inactive_means_prior_rate, digits = 3), .self$shared_variance_prior_min, .self$shared_variance_prior_max,
                                              round(.self$threshold_i, digits = 2), paste(round(.self$threshold_a, digits = 2), collapse=", "))))
    }else{
      out_prior_settigs <- data.frame(cbind(c("s0_mu", "s0_sigma", "s1_shape", "s1_rate", "tau_rate", "tau_shape", "alpha_r_shape", "alpha_r_rate",
                                              "weight_active_shape_1", "weight_active_shape_2", "weight_within_active_alpha", "spike_prior_shape_1", "spike_prior_shape_2",
                                              "active_means_dif_prior_shape", "active_meaans_dif_prior_rate", "active_variances_prior_min", "active_variances_prior_max",
                                              "inactive_means_prior_shape", "inactive_means_prior_rate", "inactive_variances_prior_min", "inactive_variances_prior_max",
                                              "threshold_i", "threshold_a"),
                                            c(.self$s0_mu, .self$s0_sigma, .self$s1_shape, .self$s1_rate, .self$tau_rate, .self$tau_shape,  .self$alpha_r_shape, .self$alpha_r_rate,
                                              .self$weight_active_shape_1, .self$weight_active_shape_2, .self$weight_within_active_alpha[1], .self$spike_prior_shape_1, .self$spike_prior_shape_2,
                                              .self$active_means_dif_prior_shape, round(.self$active_means_dif_prior_rate, digits = 3) , .self$active_variances_prior_min, .self$active_variances_prior_max,
                                              .self$inactive_means_prior_shape, round(.self$inactive_means_prior_rate, digits = 3), .self$inactive_variances_prior_min, .self$inactive_variances_prior_max,
                                              round(.self$threshold_i, digits = 2), paste(round(.self$threshold_a, digits = 2), collapse=", "))))

    }


    write.table(out_prior_settigs, file = paste0(output_directory, "/", "hyperparameter_settings.txt"),
                sep = "\t", row.names = F, quote = F, col.names = F)


    ###### YgLikelihood is not used ???? delete?
    .self$YgLikelihood[inactive_idx] <- .self$computeInactiveLikelihood(y = Yg[inactive_idx], spike_p = spike_probability, im = inactive_means, iv = inactive_variances, spike_allocation = inactive_spike_allocation[inactive_idx])
    .self$YgLikelihood[active_idx] <- .self$computeActiveLikelihood(y = Yg[active_idx], active_component = allocation_within_active[[1]][active_idx], am = active_means, av = active_variances)


    #######################
    #####Book keeping######
    #######################

    .self$field_list_values <- list(allocation_active_inactive, allocation_within_active, inactive_spike_allocation,
                                    weight_active, weight_within_active,
                                    inactive_means, inactive_variances, active_means, active_means_dif, active_variances,
                                    Xg, Yg, s0, s1, Sg, variance_g, alpha_r, p_x,
                                    active_idx, inactive_idx, in_spike_idx, out_spike_idx)

    .self$field_list_names <- list("allocation_active_inactive", "allocation_within_active", "inactive_spike_allocation",
                                   "weight_active", "weight_within_active",
                                   "inactive_means", "inactive_variances", "active_means", "active_means_dif", "active_variances",
                                   "Xg", "Yg", "s0", "s1", "Sg", "variance_g", "alpha_r", "p_x",
                                   "active_idx", "inactive_idx", "in_spike_idx", "out_spike_idx")


    .self$lnl_trace[[1]] <- sum(XgLikelihood)


  }

)
