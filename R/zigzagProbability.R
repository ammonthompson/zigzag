zigzag$methods(

  computeTauPriorProbability = function(x){

    lnp <- dgamma(x, tau_shape, tau_rate, log = TRUE)

    return(lnp)

  },

  computeVarianceGPriorProbability = function(x, sx, gtau = tau){

    # lnp <- dlnorm(x, log(sx) + gtau, sqrt(gtau), log=TRUE)
    # lnp <- dgamma(x, shape = sx/gtau, rate = 1/gtau, log = TRUE)

    rtgtau = sqrt(gtau)

    if(variance_g_upper_bound == Inf){

      lnp <- -log(x * rtgtau * sqrt2pi) - 0.5 * ((log(x) - log(sx) - gtau)/rtgtau)^2

    }else{

      lnp <- -log(x * rtgtau * sqrt2pi) - 0.5 * ((log(x) - log(sx) - gtau)/rtgtau)^2 -
        plnorm(variance_g_upper_bound, log(sx) + gtau, sqrt(gtau), log.p = TRUE)

    }

    return(lnp)

  },

  computeS0PriorProbability = function(x) {

    lnp <- dnorm(x, s0_mu, s0_sigma, log = TRUE)

    return(lnp)

  },

  computeS1PriorProbability = function(x){

    lnp <- dgamma(abs(x), s1_shape, s1_rate, log = TRUE)

    return(lnp)

  },

  computeAlphaRPriorProbability = function(x){

    lnp <- sum(dgamma(x, alpha_r_shape, alpha_r_rate, log = TRUE))

    return(lnp)

  },

  computeBiasMatrixProbability = function(i_idx, a_idx){
    lnp <- c()
    lnp[i_idx] <- computeBiasProbability(inactive_bias)
    lnp[a_idx] <- computeBiasProbability(active_bias)

    return(lnp)
  },

  computeBiasProbability = function(x){

    dir_val = x/bias_scalor + 1/num_libraries

    log_jacobian = -num_libraries * log(bias_scalor)

    lnp = log(d_dirichlet(dir_val, rep(bias_alpha, num_libraries))) + log_jacobian

    return(lnp)

  },

  computeXgLikelihood = function(xg, yg, sx, px,
                                 spike_allocation = inactive_spike_allocation,
                                 spike_prob = spike_probability,
                                 nd_spike = no_detect_spike,
                                 recover_x = FALSE,
                                 bias_matrix = lib_bias_matrix) { # returns a vecor of log-likelihoods, one for each gene x given y and Sg and p_x

    xg <- as.matrix(xg)
    px <- as.matrix(px)
    bias_matrix <- as.matrix(bias_matrix)

    inspike_idx <- which(spike_allocation == 1)
    outspike_idx <- which(spike_allocation == 0)

    lnl = NULL

    if(length(outspike_idx) > 1){

      lnl[outspike_idx] <- rowSums(as.matrix(sapply(1:ncol(xg), function(lib){

        ll <- NULL

        not_detected_idx <- which(xg[outspike_idx, lib] == inf_tol)

        detected_idx <- which(xg[outspike_idx, lib] != inf_tol)

        ll[not_detected_idx] <- log(1 - px[outspike_idx,lib][not_detected_idx])

        ll[detected_idx] <- log(px[outspike_idx,lib][detected_idx]) +
          dnorm(xg[outspike_idx,lib][detected_idx], yg[outspike_idx][detected_idx] + bias_matrix[outspike_idx,lib][detected_idx],
                sqrt(sx[outspike_idx][detected_idx]), log=TRUE)


        return(ll)

      })))

    }

    if(length(outspike_idx) == 1){

      lnl[outspike_idx] <- rowSums(t(as.matrix(sapply(1:ncol(xg), function(lib){

        ll <- NULL

        not_detected_idx <- which(xg[outspike_idx,lib] == inf_tol)

        detected_idx <- which(xg[outspike_idx,lib] != inf_tol)


        ll[not_detected_idx] <- log(1 - px[outspike_idx,lib][not_detected_idx])

        ll[detected_idx] <- log(px[outspike_idx,lib][detected_idx]) +
          dnorm(xg[outspike_idx,lib][detected_idx], yg[outspike_idx][detected_idx] + bias_matrix[outspike_idx,lib][detected_idx],
                sqrt(sx[outspike_idx][detected_idx]), log=TRUE)

        return(ll)

      }))))

    }


    lnl[inspike_idx] <- log(spike_allocation[inspike_idx] * nd_spike[inspike_idx])

    if(recover_x) recover()

    if(beta < 0.0000001){
      lnl[seq(length(lnl))] = 0
    }

    return(beta * lnl)

  },

  computeActiveWeightPriorProbability = function(x) {

    lnp <- dbeta(x, weight_active_shape_1, weight_active_shape_2, log=TRUE)

    return(lnp)

  },

  computeWithinActiveWeightPriorProbability = function(x) {
    lnp <- .self$d_dirichlet(x, weight_within_active_alpha, log = TRUE)
    return(lnp)
  },

  computeInactiveMeansPriorProbability = function(x) {

    #threshold_i - noise_expression ~ gamma
    # lnp <- sum(sapply(x, function(y) sum(dgamma((threshold_i - y), shape=inactive_means_prior_shape, rate=inactive_means_prior_rate, log=TRUE))))

    lnp <- dgamma((threshold_i - x), shape=inactive_means_prior_shape, rate=inactive_means_prior_rate, log=TRUE)

    return(lnp)

  },

  computeInactiveVariancesPriorProbability = function(x) {

    lnp <- dunif(log(x, 10), min=inactive_variances_prior_log_min, max=inactive_variances_prior_log_max, log=TRUE) #- log(x) - 0.8340324

    return(lnp)

  },

  computeSpikePriorProbability = function(x) {

    lnp <- sum(dbeta(x, shape1=spike_prior_shape_1, shape2=spike_prior_shape_2, log=TRUE))

    return(lnp)
  },

  computeActiveMeansDifPriorProbability = function(x) {

    lnp <- sum(sapply(x, function(y) sum(dgamma(y, shape=active_means_dif_prior_shape, rate=active_means_dif_prior_rate, log=TRUE))))

    return(lnp)

  },

  computeActiveVariancesPriorProbability = function(x) {

    #lnp <- sum(sapply(x, function(y) sum(dunif(log(y, 10), min=active_variances_prior_log_min, max=active_variances_prior_log_max, log=TRUE))) - log(x) - 0.8340324)
    lnp <- sum(sapply(x, function(y) sum(dunif(log(y, 10), min=active_variances_prior_log_min, max=active_variances_prior_log_max, log=TRUE))))

    return(lnp)

  },

  computeAllocationPriorProbability = function(p){ # x is the allocations {0,1,...,K} and p is the corresponding probs such that sum(p) = 1

    lnp <- sum(dbinom(length(inactive_idx), size = num_transcripts, prob = p[1], log = TRUE)) +
        sum(log(1 - p[1]) + dmultinom(tabulate(allocation_within_active[[1]][active_idx], nbins = num_active_components), prob = p[-1], log=TRUE))

    return(lnp)

  },

  computeActiveAllocationPriorProbability = function(x, p) {

    lnp <- sum(dbinom(x, 1, p, log=TRUE))

    return(lnp)

  }, #depricated

  computeWithinActiveAllocationPriorProbability = function(x, p, active) {

    lnp <- sum(dmultinom(table(x[[1]][active]), prob = p, log=TRUE))

    return(lnp)

  }, #depricated

  computeInactiveLikelihood = function(y, spike_p, im, iv, spike_allocation = inactive_spike_allocation, all = FALSE) { # returns a vecor of log-likelihoods, one for each gene x allocated to inacive component

    is <- sqrt(iv)

    if(all){

      lnl <- c()

      lnl[in_spike_idx] <- log(spike_p)

      lnl[out_spike_idx] <- log((1 - spike_p)) - log(sqrt2pi * is) - (y[out_spike_idx] - im)^2/(2 * iv)

    }else{

      inspike <- which(spike_allocation == 1)

      outspike <- which(spike_allocation == 0)

      lnl <- c()

      lnl[inspike] <- log(spike_p)

      lnl[outspike] <- log(1 - spike_p) - log(sqrt2pi * is) - (y[outspike] - im)^2/(2 * iv)

    }

    return(lnl)

  },

  computeActiveLikelihood = function(y, active_component, am, av, spike_allocation = inactive_spike_allocation, all = FALSE) { # returns a vecor of log-likelihoods, one for each gene x allocated to an active component

    if(all) active_component <- rep(active_component, length(y))

    if(is.list(active_component)){

      lnl <- NULL

      for(k in seq(num_active_components)){

        k_idx = which(active_component == k)

        lnl[k_idx] <- -log(sqrt2pi * sqrt(av[k])) - (y[k_idx] - am[k])^2/(2 * av[k])

      }


    }else{

      lnl <- NULL

      for(k in seq(num_active_components)){

        k_idx = which(active_component == k)

        lnl[k_idx] <- -log(sqrt2pi * sqrt(av[k])) - (y[k_idx] - am[k])^2/(2 * av[k])

      }

    }

    return(lnl)

  },



  ###############
  ### Develop ###
  ###############

  old_computeInactiveLikelihood = function(y, spike_p, im, iv, spike_allocation = inactive_spike_allocation) { # returns a vecor of log-likelihoods, one for each gene x allocated to inacive component

    is <- sqrt(iv)

    lnl <- beta * log(spike_p * spike_allocation  + (1 - spike_p) * (1 - spike_allocation) * dnorm(y, im, is))

    return(lnl)

  },

  old2_computeInactiveLikelihood = function(y, spike_p, im, iv, spike_allocation = inactive_spike_allocation, all = FALSE) { # returns a vecor of log-likelihoods, one for each gene x allocated to inacive component

    is <- sqrt(iv)

    if(all){

      lnl <- c()

      lnl[in_spike_idx] <- log(spike_p)

      lnl[out_spike_idx] <- (log((1 - spike_p)) + dnorm(y[out_spike_idx], im, is, log = TRUE))

    }else{

      inspike <- which(spike_allocation == 1)

      outspike <- which(spike_allocation == 0)

      lnl <- c()

      lnl[inspike] <- log(spike_p)

      lnl[outspike] <- (log((1 - spike_p)) + dnorm(y[outspike], im, is, log = TRUE))

    }

    return(lnl)

  },

  old1_computeActiveLikelihood = function(y, active_component, am, av, spike_allocation = inactive_spike_allocation) { # returns a vecor of log-likelihoods, one for each gene x allocated to an active component

    if(is.list(active_component)){

      lnl <- dnorm(y, am[ active_component[[1]] ], sqrt(av[ active_component[[1]] ]), log = TRUE)

    }else{

      lnl <- dnorm(y, am[active_component], sqrt(av[active_component]), log = TRUE)

    }

    return(lnl)

  },

  computeXgLikelihood_oneLibrary = function(xg, px,
                                            spike_allocation = inactive_spike_allocation,
                                            spike_prob = spike_probability,
                                            nd_spike = no_detect_spike, recover_x = FALSE) { # returns a vecor of log-likelihoods, one for each gene x given y and Sg and p_x

    xg <- as.matrix(xg)
    px <- as.matrix(px)

    inspike_idx <- which(spike_allocation == 1)
    outspike_idx <- which(spike_allocation == 0)

    lnl = rep(NA, length(px))

    if(length(outspike_idx) > 0){

      not_detected_idx <- which(xg[outspike_idx] == inf_tol)

      detected_idx <- which(xg[outspike_idx] != inf_tol)

      detected_inactive_idx <- which(allocation_active_inactive[outspike_idx][detected_idx] == 0)

      detected_active_idx <- which(allocation_active_inactive[outspike_idx][detected_idx] == 1)

      lnl[outspike_idx][not_detected_idx] <- log(1 - px[outspike_idx][not_detected_idx])

      if(recover_x) recover()

      if(length(detected_inactive_idx) > 0){

        # lnl[outspike_idx][detected_idx][detected_inactive_idx] <- log(px[outspike_idx][detected_idx][detected_inactive_idx]) +
        #   log((1 - weight_active) * exp(.self$computeInactiveLikelihood(Yg[outspike_idx][detected_idx][detected_inactive_idx], spike_probability,
        #                                                                 inactive_means, inactive_variances, inactive_spike_allocation[outspike_idx][detected_idx][detected_inactive_idx])))
        #

        lnl[outspike_idx][detected_idx][detected_inactive_idx] <- log(px[outspike_idx][detected_idx][detected_inactive_idx])
      }

      if(length(detected_active_idx) > 0){


        # lnl[outspike_idx][detected_idx][detected_active_idx] <- log(px[outspike_idx][detected_idx][detected_active_idx]) +
        #   log( weight_active * exp(.self$computeActiveLikelihood(Yg[outspike_idx][detected_idx][detected_active_idx],
        #                                                          allocation_within_active[[1]][outspike_idx][detected_idx][detected_active_idx], active_means, active_variances)))


        lnl[outspike_idx][detected_idx][detected_active_idx] <- log(px[outspike_idx][detected_idx][detected_active_idx])

      }

    }

    lnl[inspike_idx] <- log(spike_allocation[inspike_idx] * nd_spike[inspike_idx])

    return(beta * lnl)

  },

  test_parallel_lxg = function(xg, yg, sx, px, spike_allocation = inactive_spike_allocation, spike_prob = spike_probability, nd_spike = no_detect_spike, num_cores = detectCores()){

    cl <- makeCluster(num_cores)

    lnl <- rowSums(matrix(unlist(clusterApply(cl, 1:ncol(xg), fun = function(lib){

      return(.self$computeXgLikelihood(xg[,lib], yg, sx, px[,lib], spike_allocation, spike_prob, nd_spike))
    })), ncol = ncol(xg), byrow= FALSE))

    stopCluster(cl)

    return(lnl)


  }

)











