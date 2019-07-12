zigzag$methods(

  ##########################
  ### Lower level params ###
  ##########################

  mhSigma_g = function(recover_x = FALSE, tune = FALSE){

    HR = 0

    proposal <- .self$scale_move(sigma_g[out_spike_idx], tuningParam_sigma_g[out_spike_idx])

    sigma_g_proposed <- proposal[[1]]

    HR <- log(proposal[[2]])

    proposed_XgLikelihood <- .self$computeXgLikelihood(Xg[out_spike_idx,], Yg[out_spike_idx], sigma_g_proposed, p_x[out_spike_idx,],
                                                       spike_allocation = inactive_spike_allocation[out_spike_idx], nd_spike = no_detect_spike[out_spike_idx])

    proposed_sigmaGProbability <- .self$computeSigmaGPriorProbability(sigma_g_proposed, Sg[out_spike_idx])

    Lp_pp = proposed_XgLikelihood + proposed_sigmaGProbability

    Lc_pp = XgLikelihood[out_spike_idx] + sigma_g_probability[out_spike_idx]

    R = temperature * (Lp_pp - Lc_pp) + HR

    if(recover_x) recover()

    current_and_proposed <- cbind(sigma_g[out_spike_idx], sigma_g_proposed)

    current_and_proposed_likelihood <- cbind(XgLikelihood[out_spike_idx], proposed_XgLikelihood)


    num_choices <- length(out_spike_idx)

    choice_idx <- (log(runif(num_choices)) < R) + 1  # 1 means reject proposal, 2 means accept proposal.

    choice_matrix <- cbind(seq(num_choices), choice_idx) # row-column pairs for current_and_proposed_sigma_g_probability with 1st column meaning reject proposal

    sigma_g[out_spike_idx] <<- current_and_proposed[choice_matrix]

    if(tune){

      if(length(out_spike_idx) < 2){

        sigma_g_trace[out_spike_idx,] <<- cbind(matrix(sigma_g_trace[out_spike_idx,-1], nrow = length(out_spike_idx)), choice_idx - 1)

      }else{

        sigma_g_trace[out_spike_idx,] <<- cbind(sigma_g_trace[out_spike_idx,-1], choice_idx - 1)
      }
    }


    XgLikelihood[out_spike_idx] <<- current_and_proposed_likelihood[choice_matrix]

    current_and_proposed_sigma_g_probability <- cbind(sigma_g_probability[out_spike_idx], proposed_sigmaGProbability)

    sigma_g_probability[out_spike_idx] <<- current_and_proposed_sigma_g_probability[choice_matrix]

  },

  mhSg = function(recover_x = FALSE, tune = FALSE){

    HR = 0

    selection = sample(seq(2), 1, prob = c(1,1))


    if( selection == 1){

      proposed_s0 <- s0 + rnorm(1, 0, tuningParam_s0)

      proposed_s1 <- s1

      if(tune) s0_trace[[1]][[1]] <<- c(s0_trace[[1]][[1]][-1], 0)

    }else{

      proposed_s0 <- s0

      proposal <- .self$scale_move(s1, tuningParam_s1)

      proposed_s1 <- proposal[[1]]

      HR <- log(proposal[[2]])

      if(tune) s1_trace[[1]][[1]] <<- c(s1_trace[[1]][[1]][-1], 0)

    }

    proposed_Sg <- .self$get_sigmax(ss0 = proposed_s0, ss1 = proposed_s1)

    proposed_sigma_g_Likelihood <- .self$computeSigmaGPriorProbability(sigma_g[out_spike_idx], proposed_Sg[out_spike_idx])

    Lp <- sum(proposed_sigma_g_Likelihood)

    Lc <- sum(sigma_g_probability[out_spike_idx])

    pp <- .self$computeS0PriorProbability(proposed_s0) + .self$computeS1PriorProbability(proposed_s1)

    pc <- .self$computeS0PriorProbability(s0) + .self$computeS1PriorProbability(s1)

    R <- temperature * (Lp - Lc + pp - pc) + HR

    if(recover_x) recover()


    if(log(runif(1)) < R){

      s0 <<- proposed_s0
      s1 <<- proposed_s1

      if(selection == 1){

        if(tune) s0_trace[[1]][[1]][100] <<- 1

      }else {

        if(tune) s1_trace[[1]][[1]][100] <<- 1

      }

      Sg <<- proposed_Sg

      sigma_g_probability[out_spike_idx] <<- proposed_sigma_g_Likelihood

    }

  },

  mhTau = function(recover_x = FALSE, tune = FALSE){

    # proposed_tau <- abs(tau + rnorm(1, 0, tuningParam_tau))

    proposal <- .self$scale_move(tau, tuningParam_tau)

    proposed_tau <- proposal[[1]]

    HR <- log(proposal[[2]])

    proposed_tau_Likelihood <- .self$computeSigmaGPriorProbability(sigma_g[out_spike_idx], Sg[out_spike_idx], gtau = proposed_tau)

    # tau_Likelihood <- .self$computeSigmaGPriorProbability(sigma_g[out_spike_idx], Sg[out_spike_idx], gtau = tau)

    tau_Likelihood <- sigma_g_probability[out_spike_idx]

    Lp <- sum(proposed_tau_Likelihood)
    Lc <- sum(tau_Likelihood)

    pp <- .self$computeTauPriorProbability(proposed_tau)
    pc <- .self$computeTauPriorProbability(tau)

    R <- temperature * (Lp - Lc + pp - pc) + HR

    if(recover_x) recover()

    if(tune) tau_trace[[1]][[1]] <<- c(tau_trace[[1]][[1]][-1], 0)

    if(log(runif(1)) < R){

      tau <<- proposed_tau

      if(tune) tau_trace[[1]][[1]][100] <<- 1

      sigma_g_probability[out_spike_idx] <<- proposed_tau_Likelihood

    }


  },

  mhS0Tau = function(recover_x = FALSE, tune = FALSE){
    HR = 0

    proposal <- .self$scale_move(tau, tuningParam_s0tau)

    proposed_s0 <- s0 - log(proposal[[2]])

    proposed_tau <- proposal[[1]]

    HR <- log(proposal[[2]])

    if(tune) s0tau_trace[[1]][[1]] <<- c(s0tau_trace[[1]][[1]][-1], 0)

    proposed_Sg <- .self$get_sigmax(ss0 = proposed_s0, ss1 = s1)

    proposed_sigma_g_Likelihood <- .self$computeSigmaGPriorProbability(sigma_g[out_spike_idx], proposed_Sg[out_spike_idx], proposed_tau)

    Lp <- sum(proposed_sigma_g_Likelihood)

    Lc <- sum(sigma_g_probability[out_spike_idx])

    pp <- .self$computeS0PriorProbability(proposed_s0) + .self$computeTauPriorProbability(proposed_tau)

    pc <- .self$computeS0PriorProbability(s0) + .self$computeTauPriorProbability(tau)

    R <- temperature * (Lp - Lc + pp - pc) + HR

    if(recover_x) recover()

    if(log(runif(1)) < R){

      s0 <<- proposed_s0

      tau <<- proposed_tau

      Sg <<- proposed_Sg

      if(tune) s0tau_trace[[1]][[1]][100] <<- 1

      sigma_g_probability[out_spike_idx] <<- proposed_sigma_g_Likelihood

    }

  },

  mhP_x = function(recover_x = FALSE, tune = FALSE){

    randlib=sample(num_libraries,1)

    proposed_alpha_r <- alpha_r

    proposal <- .self$scale_move(alpha_r[randlib], tuningParam_alpha_r[randlib])

    proposed_alpha_r[randlib] <- proposal[[1]]

    if(tune) alpha_r_trace[[1]][[1]][randlib,] <<- c(alpha_r_trace[[1]][[1]][randlib, -1], 0)

    HR <- log(proposal[[2]])

    proposed_p_x <- p_x

    if(recover_x) recover()

    proposed_p_x[,randlib] <- .self$get_px(falpha_r = proposed_alpha_r[randlib], yy = Yg)
    # proposed_p_x[,randlib] <- .self$xxget_px(falpha_r = proposed_alpha_r[randlib], xx = Xg[,randlib])

    if(num_libraries > 2){

      proposed_xgLikelihood <- XgLikelihood + .self$computeXgLikelihood(Xg[,randlib], Yg, sigma_g, proposed_p_x[,randlib]) -
        .self$computeXgLikelihood(Xg[,randlib], Yg, sigma_g, p_x[,randlib])

    }else{

      proposed_xgLikelihood <- .self$computeXgLikelihood(Xg, Yg, sigma_g, proposed_p_x)

    }

    Lp <- sum(proposed_xgLikelihood)

    Lc <- sum(XgLikelihood)

    pp <- .self$computeAlphaRPriorProbability(proposed_alpha_r)

    pc <- .self$computeAlphaRPriorProbability(alpha_r)

    R <- temperature * (Lp - Lc + pp - pc) + HR


    if(log(runif(1)) < R & R != Inf){

      alpha_r[randlib] <<- proposed_alpha_r[randlib]

      p_x[,randlib] <<- proposed_p_x[,randlib]

      if(tune) alpha_r_trace[[1]][[1]][randlib,100] <<- 1

      XgLikelihood <<- proposed_xgLikelihood

    }

  },


  ##############################
  ### Upper level parameters ###
  ##############################

  mhYg = function(recover_x = FALSE, tune = FALSE){

    Yg_proposed <<- Yg + rnorm(num_transcripts, 0, tuningParam_yg) * (1 - inactive_spike_allocation) #only change out of spike genes

    Sg_proposed <- .self$get_sigmax(yy = Yg_proposed)

    p_x_proposed <- .self$get_px(yy = Yg_proposed)

    ## proposed Yg likelihoods
    proposed_XgLikelihood = .self$computeXgLikelihood(Xg, Yg_proposed, sigma_g, p_x_proposed)

    proposed_sigmaGProbability = .self$computeSigmaGPriorProbability(sigma_g, Sg_proposed)

    Lp_pp = proposed_XgLikelihood

    Lp_pp[inactive_idx] = Lp_pp[inactive_idx] +
      .self$computeInactiveLikelihood(Yg_proposed[inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx])

    Lp_pp[active_idx] = Lp_pp[active_idx] +
      .self$computeActiveLikelihood(Yg_proposed[active_idx],  allocation_within_active[[1]][active_idx], active_means, active_variances)

    Lp_pp = Lp_pp + proposed_sigmaGProbability


    ## current Yg likelihoods
    Lc_pc = XgLikelihood

    Lc_pc[inactive_idx] = Lc_pc[inactive_idx] +
      .self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx])

    Lc_pc[active_idx] = Lc_pc[active_idx] +
      .self$computeActiveLikelihood(Yg[active_idx],  allocation_within_active[[1]][active_idx], active_means, active_variances)

    Lc_pc = Lc_pc + sigma_g_probability

    R = temperature * (Lp_pp - Lc_pc)

    current_and_proposed <- cbind(Yg, Yg_proposed)

    current_and_proposed_likelihood <- cbind(XgLikelihood, proposed_XgLikelihood)

    if(recover_x) recover()

    choice_idx <- 1 + (log(runif(num_transcripts)) < R)

    choice_matrix <- cbind(seq(num_transcripts), choice_idx)

    Yg <<- current_and_proposed[cbind(seq(num_transcripts), choice_idx)]

    if(tune){

      if(length(out_spike_idx) < 2){

        Yg_trace[out_spike_idx,] <<- cbind(matrix(Yg_trace[out_spike_idx,-1], nrow = length(out_spike_idx)), choice_idx[out_spike_idx] - 1)

      }else{

        Yg_trace[out_spike_idx,] <<- cbind(Yg_trace[out_spike_idx,-1], choice_idx[out_spike_idx] - 1)
      }
    }

    ## update XgLikelihood, Sg and p_x

    XgLikelihood <<- current_and_proposed_likelihood[choice_matrix]

    .self$set_sigmaX_pX()

    current_and_proposed_sigmaGProbability <- cbind(sigma_g_probability, proposed_sigmaGProbability)

    sigma_g_probability[out_spike_idx] <<- current_and_proposed_sigmaGProbability[choice_matrix][out_spike_idx]

  },

  gibbsMixtureWeights = function(recover_x = FALSE, tune = FALSE){
    #Propose new weights for the vector: c(1 - omega^a, omega^a * c(omega^m_1/omega_a, omega^m_2, ...))

    all_allocations = allocation_within_active[[1]] * allocation_active_inactive

    num_in_components <- sapply(0:num_active_components, function(y){return(sum(all_allocations == y))}) * beta

    new_weights <- r_dirichlet(c(weight_active_shape_2, weight_within_active_alpha) + num_in_components)

    if(recover_x) recover()

    if(! is.nan(new_weights[1])){

      weight_active <<- 1 - new_weights[1]

      weight_within_active <<- new_weights[-1]/sum(new_weights[-1]) #conditioned on omega^a (sum of the new_weights minus that of the inactive component)

    }

  },

  gibbsAllocationActiveInactive = function(recover_x = FALSE, tune = FALSE) {

    ### Gibbs move active/inactive component

    L_0 = exp(.self$computeInactiveLikelihood(Yg, spike_probability, inactive_means, inactive_variances, inactive_spike_allocation, all = TRUE))

    L_1 = 0

    for(k in seq(num_active_components)){

      L_1 = L_1 + weight_within_active[k] * exp(.self$computeActiveLikelihood(Yg, k, active_means, active_variances, all = TRUE))

    }

    prob_active <- weight_active * L_1 * (1-inactive_spike_allocation)/((1 - weight_active) * L_0 + weight_active * L_1 * (1-inactive_spike_allocation))

    prob_active[is.nan(prob_active)] <- weight_active
    allocation_active_inactive <<- rbinom(num_transcripts, 1, prob_active)

    if(!is.null(active_gene_set)) allocation_active_inactive[active_gene_set_idx] <<- as.integer(1)

    .self$setActiveInactive_idx()

  },

  gibbsAllocationWithinActive = function(recover_x = FALSE, tune = FALSE){

    ln_smallL <- 0 #-400 #probably delete this

    if(length(active_idx) > 0){

      logL_k = NULL

      if(length(active_idx) > 1){

        for(k in seq(num_active_components)){

          logL_k = cbind(logL_k, log(weight_within_active[k]) +
                           .self$computeActiveLikelihood(Yg[active_idx], k, active_means, active_variances, all = TRUE) - ln_smallL)

        }

      }else{

        for(k in seq(num_active_components)){

          logL_k = c(logL_k, log(weight_within_active[k]) +
                       .self$computeActiveLikelihood(Yg[active_idx], k, active_means, active_variances, all = TRUE) - ln_smallL)

        }

      }

      if(is.vector(logL_k)) logL_k = matrix(logL_k, nrow=1, ncol=num_active_components)

      if(recover_x) recover()

      L_k = exp(logL_k)
      # L_k[which(is.na(L_k), arr.ind=T)]=0

      probs = L_k/rowSums(L_k)
      cumprobs = rowCumsums(probs)
      randsamp = floor(cumprobs/runif(nrow(probs)))
      zeros = matrix(0, nrow(probs), num_active_components)
      new_allocation <- rowSums(randsamp==zeros) + 1
      newallo <- which(! is.na(new_allocation))
      allocation_within_active[[1]][active_idx][newallo] <<- new_allocation[newallo]

    }



  },

  mhInactiveMeans = function(recover_x = FALSE, tune = FALSE){
    ### for each distribution parameter: propose and accept new value with MH algo

    inactive_means_proposed <<- threshold_i - abs(threshold_i - inactive_means - rnorm(1, 0, inactive_mean_tuningParam[[1]]))

    pp = .self$computeInactiveMeansPriorProbability(inactive_means_proposed)

    pc = .self$computeInactiveMeansPriorProbability(inactive_means)

    if(length(inactive_idx) > 0){

      Lp = sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means_proposed, inactive_variances,
                                               inactive_spike_allocation[inactive_idx]))
      Lc = sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances,
                                               inactive_spike_allocation[inactive_idx]))

    }else{

      Lp = 0
      Lc = 0

    }

    R = temperature * (Lp - Lc + pp - pc)


    if(tune) inactive_means_trace[[1]][[1]] <<- c(inactive_means_trace[[1]][[1]][-1],0)

    if(recover_x) recover()

    if(log(runif(1)) < R){

      inactive_means <<- inactive_means_proposed

      if(tune) inactive_means_trace[[1]][[1]][100] <<- 1

    }

  },

  mhInactiveVariances = function(recover_x = FALSE, tune = FALSE){

    inactive_variances_proposed <<- 10^(.self$two_boundary_slide_move(log(inactive_variances, 10), inactive_variances_prior_log_min,
                                                                      inactive_variances_prior_log_max, inactive_variance_tuningParam))
    if(tune) inactive_variances_trace[[1]][[1]] <<- c(inactive_variances_trace[[1]][[1]][-1], 0)

    pp = .self$computeInactiveVariancesPriorProbability(inactive_variances_proposed)

    pc = .self$computeInactiveVariancesPriorProbability(inactive_variances)

    if(length(inactive_idx) > 0){

      Lp = sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances_proposed,
                                               inactive_spike_allocation[inactive_idx]))
      Lc = sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances,
                                               inactive_spike_allocation[inactive_idx]))

    }else{

      Lp = 0
      Lc = 0

    }

    R =  temperature * (Lp - Lc + pp - pc)

    if(recover_x) recover()

    if(log(runif(1)) < R){

      inactive_variances <<- inactive_variances_proposed

      if(tune) inactive_variances_trace[[1]][[1]][100] <<- 1

    }

  },

  gibbsSpikeProb = function(recover_x = FALSE, tune = FALSE){

    num_inspike = length(in_spike_idx)
    num_outspike = length(inactive_idx) - num_inspike

    new_spike_probability <- rbeta(1, num_inspike + spike_prior_shape_1, num_outspike + spike_prior_shape_2)

    if(recover_x) recover()

    if(! is.nan(new_spike_probability)){

      spike_probability <<- new_spike_probability

    }

  },

  mhSpikeAllocation = function(recover_x = FALSE, tune = FALSE){

    # note, proposed moves are only within the inactive allocation,
    # no moves between the spike and the active component are ever proposed

    zero_and_inactive_idx = all_zero_idx[allocation_active_inactive[all_zero_idx] == 0] #this is to speed things up (will mess with simulating under prior)

    inactive_spike_allocation_proposal = abs(inactive_spike_allocation[zero_and_inactive_idx] - 1)
    proposed_out_spike_idx = which(inactive_spike_allocation_proposal==0)
    proposed_in_spike_idx = which(inactive_spike_allocation_proposal==1)

    current_out_spike_idx = which(inactive_spike_allocation[zero_and_inactive_idx]==0)
    current_in_spike_idx = which(inactive_spike_allocation[zero_and_inactive_idx]==1)


    #*** Make proposals. If in spike then propose out-of-spike inactive expression value. If out of spike propose into spike

    yg_proposed = rnorm(length(zero_and_inactive_idx), inactive_means, sqrt(inactive_variances))

    Sg_proposed <- .self$get_sigmax(yy = yg_proposed)

    sigma_g_proposed <- rlnorm(length(zero_and_inactive_idx), log(Sg_proposed) + tau, sqrt(tau))  # P(sigma_g | proposed_Yg)

    p_x_proposed <- .self$get_px(yy = yg_proposed, gl = gene_lengths[zero_and_inactive_idx])


    #*** Compute prior probabilities

    p_sigma_g_proposed <- .self$computeSigmaGPriorProbability(sigma_g_proposed, Sg_proposed)

    p_sigma_g <- .self$computeSigmaGPriorProbability(sigma_g[zero_and_inactive_idx], Sg[zero_and_inactive_idx])

    p_yg_proposed <- .self$computeInactiveLikelihood(Yg_proposed, spike_probability, inactive_means,
                                                     inactive_variances, inactive_spike_allocation_proposal)

    p_yg <- .self$computeInactiveLikelihood(Yg[zero_and_inactive_idx], spike_probability, inactive_means,
                                            inactive_variances, inactive_spike_allocation[zero_and_inactive_idx])

    pp = NULL; pc = NULL

    ### this needs a clear explanation.
    pp[c(proposed_out_spike_idx, proposed_in_spike_idx)] = c(p_yg_proposed[proposed_out_spike_idx] + p_sigma_g_proposed[proposed_out_spike_idx],
                                                             rep(log(spike_probability), length(proposed_in_spike_idx)))

    pc[c(current_in_spike_idx, current_out_spike_idx)] = c(rep(log(spike_probability), length(current_in_spike_idx)),
                                                           p_yg[current_out_spike_idx] + p_sigma_g[current_out_spike_idx])



    #*** Compute Likelihoods

    Lp = .self$computeXgLikelihood(Xg[zero_and_inactive_idx,], yg_proposed, sigma_g_proposed, p_x_proposed,
                                   inactive_spike_allocation_proposal, nd_spike = no_detect_spike[zero_and_inactive_idx])
    Lc = XgLikelihood[zero_and_inactive_idx]

    both_inf_boolean = ifelse(Lp == -Inf & Lc == -Inf, TRUE, FALSE)
    Lp[both_inf_boolean] = 0
    Lc[both_inf_boolean] = 0

    #/// Compute Hastings Ratio

    HR = NULL

    # subtract out spike_probability because of use in computeInactiveLikelihood
    HR[c(proposed_out_spike_idx, proposed_in_spike_idx)] = c(-(p_yg_proposed[proposed_out_spike_idx] - log(1-spike_probability) +
                                                                 p_sigma_g_proposed[proposed_out_spike_idx]),
                                                             p_yg[proposed_in_spike_idx] - log(1-spike_probability) +
                                                               p_sigma_g[proposed_in_spike_idx])



    #*** Update Parameters

    R = temperature * (Lp - Lc + pp - pc) + HR

    choice <- (log(runif(length(zero_and_inactive_idx))) < R) + 1

    choice[which(is.na(choice) | is.nan(choice))] <- 1 #Do not update (set choice to current value) if NA or NaN.

    choice_matrix <- cbind(seq(length(zero_and_inactive_idx)), choice)

    current_proposed = cbind(inactive_spike_allocation[zero_and_inactive_idx], inactive_spike_allocation_proposal)

    if(recover_x) recover()


    inactive_spike_allocation[zero_and_inactive_idx] <<- current_proposed[choice_matrix]

    yg_current_proposed = cbind(Yg[zero_and_inactive_idx], yg_proposed)

    Yg[zero_and_inactive_idx][proposed_out_spike_idx] <<- yg_current_proposed[choice_matrix][proposed_out_spike_idx]

    sigma_g_current_proposed <- cbind(sigma_g[zero_and_inactive_idx], sigma_g_proposed)

    sigma_g[zero_and_inactive_idx][proposed_out_spike_idx] <<- sigma_g_current_proposed[choice_matrix][proposed_out_spike_idx]

    .self$set_sigmaX_pX()

    .self$setInSpike_idx()

    xgl_current_proposed = cbind(Lc, Lp)

    XgLikelihood[zero_and_inactive_idx] <<- xgl_current_proposed[choice_matrix]

    sigma_gProb_current_proposed <- cbind(p_sigma_g, p_sigma_g_proposed)

    sigma_g_probability[zero_and_inactive_idx][proposed_out_spike_idx] <<- sigma_gProb_current_proposed[choice_matrix][proposed_out_spike_idx]

  },

  mhActiveMeansDif = function(recover_x = FALSE, tune = FALSE){
    # active means

    if(length(active_idx) > 0){

      sapply(sample(num_active_components, num_active_components, replace = TRUE), function(k){

        active_means_dif_proposed <<- active_means_dif
        active_means_dif_proposed[k] <<- abs(active_means_dif[k] + rnorm(1,0,active_mean_tuningParam[k]))

        mp = .self$calculate_active_means(active_means_dif_proposed)
        mc = active_means

        Lp <- sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], mp, active_variances, inactive_spike_allocation[active_idx]))
        Lc <- sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], mc, active_variances, inactive_spike_allocation[active_idx]))

        pp = computeActiveMeansDifPriorProbability(active_means_dif_proposed[k])
        pc = computeActiveMeansDifPriorProbability(active_means_dif[k])

        R = temperature * (Lp + pp - Lc - pc)

        if(recover_x) recover()


        if(tune) active_means_trace[[1]][[1]][k,] <<- c(active_means_trace[[1]][[1]][k,-1],0)

        if(log(runif(1)) < R){

          active_means_dif[k] <<- active_means_dif_proposed[k]
          if(tune) active_means_trace[[1]][[1]][k,100] <<- 1
          active_means <<- .self$calculate_active_means(active_means_dif_proposed)

        }

      })

    }else{


      sapply(sample(num_active_components, num_active_components, replace = TRUE), function(k){

        active_means_dif_proposed <<- active_means_dif
        active_means_dif_proposed[k] <<- abs(active_means_dif[k] + rnorm(1,0,active_mean_tuningParam[k]))

        pp = computeActiveMeansDifPriorProbability(active_means_dif_proposed[k])
        pc = computeActiveMeansDifPriorProbability(active_means_dif[k])

        R = temperature * (pp - pc)

        if(recover_x) recover()

        if(tune) active_means_trace[[1]][[1]][k,] <<- c(active_means_trace[[1]][[1]][k,-1],0)

        if(log(runif(1)) < R){

          active_means_dif[k] <<- active_means_dif_proposed[k]
          if(tune) active_means_trace[[1]][[1]][k,100] <<- 1
          active_means <<- .self$calculate_active_means(active_means_dif_proposed)

        }

      })

    }

  },

  mhActiveVariances = function(recover_x = FALSE, tune = FALSE){

    k = sample(num_active_components, num_active_components, replace = TRUE)

    if(shared_active_variances){

      active_variances_proposed <<- rep(10^(.self$two_boundary_slide_move(
        log(active_variances[1], 10), active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[1])),
        num_active_components)

      if(tune) active_variances_trace[[1]][[1]][1,] <<- c(active_variances_trace[[1]][[1]][1,-1],0)

    }else{

      active_variances_proposed[k] <<- 10^(.self$two_boundary_slide_move(log(active_variances[k], 10),
                                                                         active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[k]))
      if(tune) active_variances_trace[[1]][[1]][k,] <<- c(active_variances_trace[[1]][[1]][k,-1],0)

    }

    if(length(active_idx) > 0){

      Lp = sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], active_means, active_variances_proposed, inactive_spike_allocation[active_idx]))
      Lc = sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], active_means, active_variances, inactive_spike_allocation[active_idx]))

    }else{

      Lp = 0
      Lc = 0

    }

    pp = .self$computeActiveVariancesPriorProbability(active_variances_proposed[k])
    pc = .self$computeActiveVariancesPriorProbability(active_variances[k])

    R = temperature * (Lp + pp - Lc - pc )

    if(recover_x) recover()


    if(log(runif(1)) < R){

      if(shared_active_variances){

        active_variances <<- rep(active_variances_proposed[1], num_active_components)

        if(tune) active_variances_trace[[1]][[1]][1,100] <<- 1


      }else{

        active_variances[k] <<- active_variances_proposed[k]

        if(tune) active_variances_trace[[1]][[1]][k,100] <<- 1

      }

    }

  },




  ################################################
  # extra mcmc moves and moves Under Development #
  ################################################

  mhSpikeProb = function(recover_x = FALSE, tune = FALSE){
    # inactive spike probs

    if(length(inactive_idx) > 0){

      # spike_probability_proposed <<- abs(spike_probability + rnorm(1, 0, spike_probability_tuningParam[[1]]))
      # if(spike_probability_proposed > 1) spike_probability_proposed <<- 2 - spike_probability_proposed

      proposal_scale <- spike_probability_tuningParam

      spike_probability_proposed <<- rbeta(1, spike_probability * proposal_scale, (1 - spike_probability) * proposal_scale)

      HR = dbeta(spike_probability, spike_probability_proposed * proposal_scale, (1 - spike_probability_proposed) * proposal_scale, log = TRUE) -
        dbeta(spike_probability_proposed, spike_probability * proposal_scale, (1 - spike_probability) * proposal_scale, log = TRUE)

      Lp = sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability_proposed,
                                               inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx]))

      pp = .self$computeSpikePriorProbability(spike_probability_proposed)

      Lc = sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability,
                                               inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx]))

      pc = .self$computeSpikePriorProbability(spike_probability)

      R = temperature * (Lp - Lc + pp - pc) + HR

      if(recover_x) recover()

      if(tune) spike_probability_trace[[1]][[1]] <<- c(spike_probability_trace[[1]][[1]][-1], 0)

      if(log(runif(1)) < R){

        spike_probability <<- spike_probability_proposed

        if(tune) spike_probability_trace[[1]][[1]][100] <<- 1

      }


    }else{

      # spike_probability_proposed <<- abs(spike_probability + rnorm(1, 0, spike_probability_tuningParam[[1]]))
      # if(spike_probability_proposed > 1) spike_probability_proposed <<- 2 - spike_probability_proposed

      proposal_scale <- spike_probability_tuningParam

      spike_probability_proposed <<- rbeta(1, spike_probability * proposal_scale, (1 - spike_probability) * proposal_scale)

      HR = dbeta(spike_probability, spike_probability_proposed * proposal_scale, (1 - spike_probability_proposed) * proposal_scale, log = TRUE) -
        dbeta(spike_probability_proposed, spike_probability * proposal_scale, (1 - spike_probability) * proposal_scale, log = TRUE)

      R = temperature * (.self$computeSpikePriorProbability(spike_probability_proposed) -
                           .self$computeSpikePriorProbability(spike_probability)) + HR

      if(recover_x) recover()

      if(tune) spike_probability_trace[[1]][[1]] <<- c(spike_probability_trace[[1]][[1]][-1], 0)

      if(log(runif(1)) < R){

        spike_probability <<- spike_probability_proposed
        if(tune) spike_probability_trace[[1]][[1]][100] <<- 1

      }

    }

  },

  mhCombinedMu = function(recover_x = FALSE, tune = FALSE){

    if(recover_x) recover()

    if(rbinom(1,1,0.5) == 0){

      .self$mhInactiveMeans(recover_x, tune)

    }else{

      .self$mhActiveMeansDif(recover_x, tune)

    }

  },

  mhCombinedSigma = function(recover_x = FALSE, tune = FALSE){


    if(rbinom(1,1,0.5) == 0){

      .self$mhInactiveVariances(recover_x, tune)

    }else{

      .self$mhActiveVariances(recover_x, tune)

    }

  },

  mhCombinedMuAllocation = function(recover_x = FALSE, tune = FALSE){

    # propose new component parameters
    # propose/accept new allocations given the new parameters
    # calculate R to accept/reject the new component params and new allocations

    lib <- 1

    active_means_dif_proposed <<- active_means_dif

    inactive_means_proposed <<- inactive_means

    not_inspike <- which(inactive_spike_allocation == 0)

    # inspike <- which(inactive_spike_all == 1)

    active_or_inactive = rbinom(1, 1, 1/(num_active_components + 1))

    if(active_or_inactive == 1){

      inactive_means_proposed <<- threshold_i - abs(threshold_i - inactive_means - rnorm(1, 0, inactive_mean_tuningParam[[1]]))

      if(tune) inactive_means_trace[[1]][[lib]] <<- c(inactive_means_trace[[1]][[lib]][-1],0)

      pp = computeInactiveMeansPriorProbability(inactive_means_proposed)
      pc = computeInactiveMeansPriorProbability(inactive_means)


    }else{

      pk <- sample(seq(num_active_components), 1)

      active_means_dif_proposed[pk] <<- abs(active_means_dif[pk] + rnorm(1,0,active_mean_tuningParam[pk]))

      if(tune) active_means_trace[[1]][[lib]][pk,] <<- c(active_means_trace[[1]][[lib]][pk,-1],0)

      pp = computeActiveMeansDifPriorProbability(active_means_dif_proposed)
      pc = computeActiveMeansDifPriorProbability(active_means_dif)


    }

    # Lp = .self$computeInactiveLikelihood(Yg[not_inspike], spike_probability, inactive_means_proposed, inactive_variances,
    #                                      inactive_spike_allocation[not_inspike], all = FALSE)
    #
    # Lc = .self$computeInactiveLikelihood(Yg[not_inspike], spike_probability, inactive_means, inactive_variances,
    #                                      inactive_spike_allocation[not_inspike], all = FALSE)


    Lp = (1 - spike_probability) *
      1/(sqrt(inactive_variances) * sqrt2pi * exp((Yg[not_inspike] - inactive_means_proposed)^2/(2 * inactive_variances)))

    Lc = (1 - spike_probability) *
      1/(sqrt(inactive_variances) * sqrt2pi * exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances)))

    if(!is.null(active_gene_set)){

      actg_idx = which(gene_names[out_spike_idx] %in% active_gene_set)

      Lp[actg_idx] = 0

      Lc[actg_idx] = 0

    }


    mp = .self$calculate_active_means(active_means_dif_proposed)
    mc = active_means

    Lpk = 0
    Lck = 0

    if(length(active_idx) > 0){

      for(k in 1:num_active_components){

        Lpk = Lpk + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - mp[k])^2/(2 * active_variances[k])))

        Lck = Lck + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - mc[k])^2/(2 * active_variances[k])))

      }

    }

    lnLp <- beta_oneLib * sum(log((1 - weight_active) * Lp + weight_active * Lpk))
    lnLc <- beta_oneLib * sum(log((1 - weight_active) * Lc + weight_active * Lck))

    R = temperature * (lnLp + pp - lnLc - pc)

    if(recover_x) recover()


    if(log(runif(1)) < R & ! is.nan(R) & ! is.na(R)){

      # print("#______________________________ OLD _________________________________#")


      if(active_or_inactive == 1){

        inactive_means <<- inactive_means_proposed

        if(tune) inactive_means_trace[[1]][[lib]][100] <<- 1

      }else{

        active_means_dif[pk] <<- active_means_dif_proposed[pk]

        active_means <<- mp

        if(tune) active_means_trace[[1]][[lib]][pk,100] <<- 1


      }

      .self$gibbsAllocationActiveInactive()
      .self$gibbsAllocationWithinActive()

    }

  },

  mhCombinedSigmaAllocation = function(recover_x = FALSE, tune = FALSE){

    lib <- 1

    inactive_variances_proposed <<- inactive_variances

    active_variances_proposed <<- active_variances

    not_inspike <- which(inactive_spike_allocation == 0)

    active_or_inactive = rbinom(1,1,1/(num_active_components + 1))

    #HR = 0

    if(active_or_inactive == 1){

      inactive_variances_proposed <<- 10^(.self$two_boundary_slide_move(
        log(inactive_variances, 10), inactive_variances_prior_log_min, inactive_variances_prior_log_max, inactive_variance_tuningParam))

      pp = computeInactiveVariancesPriorProbability(inactive_variances_proposed)

      pc = computeInactiveVariancesPriorProbability(inactive_variances)

      if(tune) inactive_variances_trace[[1]][[lib]] <<- c(inactive_variances_trace[[1]][[lib]][-1],0)

      #HR = log(inactive_variances_proposed/inactive_variances)


    }else{

      if(shared_active_variances){

        active_variances_proposed <<- rep(10^(.self$two_boundary_slide_move(
          log(active_variances[1], 10), active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[1])),
          num_active_components)

        pp = computeActiveVariancesPriorProbability(active_variances_proposed[1])

        pc = computeActiveVariancesPriorProbability(active_variances[1])

        if(tune) active_variances_trace[[1]][[lib]][1,] <<- c(active_variances_trace[[1]][[lib]][1,-1],0)

        #HR = log(active_variances_proposed[1]/active_variances[1])


      }else{

        pk <- sample(seq(num_active_components), 1)

        active_variances_proposed[pk] <<- 10^(.self$two_boundary_slide_move(
          log(active_variances[pk], 10), active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[pk]))

        pp = computeActiveVariancesPriorProbability(active_variances_proposed)

        pc = computeActiveVariancesPriorProbability(active_variances)

        if(tune) active_variances_trace[[1]][[lib]][pk,] <<- c(active_variances_trace[[1]][[lib]][pk,-1],0)

        #HR = log(active_variances_proposed[pk]/active_variances[pk])

      }

    }


    Lp = (1 - spike_probability) * 1/(sqrt(inactive_variances_proposed) * sqrt2pi *
                                        exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances_proposed)))

    Lc = (1 - spike_probability) * 1/(sqrt(inactive_variances) * sqrt2pi *
                                        exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances)))

    if(!is.null(active_gene_set)){

      actg_idx = which(gene_names[out_spike_idx] %in% active_gene_set)

      Lp[actg_idx] = 0

      Lc[actg_idx] = 0

    }

    Lpk = 0
    Lck = 0

    if(length(active_idx) > 0){

      for(k in 1:num_active_components){

        Lpk = Lpk + weight_within_active[k] *
          1/(sqrt(active_variances_proposed[k]) * sqrt2pi * exp((Yg[not_inspike] - active_means[k])^2/(2 * active_variances_proposed[k])))

        Lck = Lck + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - active_means[k])^2/(2 * active_variances[k])))

      }

    }

    lnLp <- beta_oneLib * sum(log((1 - weight_active) * Lp + weight_active * Lpk))
    lnLc <- beta_oneLib * sum(log((1 - weight_active) * Lc + weight_active * Lck))

    R = temperature * (lnLp + pp - lnLc - pc) #+ HR

    if(recover_x) recover()

    if(log(runif(1)) < R & ! is.nan(R) & ! is.na(R)){

      # print("#______________________________ OLD _________________________________#")

      if(active_or_inactive == 1){

        inactive_variances <<- inactive_variances_proposed

        if(tune) inactive_variances_trace[[1]][[lib]][100] <<- 1

      }else{

        if(shared_active_variances){

          active_variances <<- rep(active_variances_proposed[1], num_active_components)

          if(tune) active_variances_trace[[1]][[lib]][1,100] <<- 1


        }else{

          active_variances[pk] <<- active_variances_proposed[pk]

          if(tune) active_variances_trace[[1]][[lib]][pk,100] <<- 1

        }


      }

      .self$gibbsAllocationActiveInactive()
      .self$gibbsAllocationWithinActive()

    }

  },

  mhMixtureWeights = function(recover_x = FALSE, tune = FALSE){
    ## active genes are independently allocated to active subcomponents across libraries
    dirichlet_proposal_scale <- mixture_weight_tuningParam

    mixtureWeights <- c(1 - weight_active, weight_active * weight_within_active)
    mixtureWeights_proposed <- r_dirichlet(mixtureWeights * dirichlet_proposal_scale)

    if(0 %in% mixtureWeights_proposed) break

    weight_active_proposed <<- 1 - mixtureWeights_proposed[1]
    weight_within_active_proposed <<- mixtureWeights_proposed[-1]/sum(mixtureWeights_proposed[-1])

    allweights_priorProbabilityProposed <- d_dirichlet(mixtureWeights_proposed, c(weight_active_shape_1, weight_within_active_alpha), log = TRUE)
    allweights_priorProbability <- d_dirichlet(mixtureWeights, c(weight_active_shape_2, weight_within_active_alpha) , log = TRUE) #weight_active_shape_2 is for the inactive comp.

    logHastings = d_dirichlet(mixtureWeights, mixtureWeights_proposed * dirichlet_proposal_scale , log = TRUE)  -
      d_dirichlet(mixtureWeights_proposed, mixtureWeights * dirichlet_proposal_scale , log =TRUE)

    logR = temperature * (.self$computeAllocationPriorProbability(c(1 - weight_active_proposed, weight_within_active_proposed)) + allweights_priorProbabilityProposed -
                            .self$computeAllocationPriorProbability(c(1 - weight_active, weight_within_active)) - allweights_priorProbability) + logHastings


    if(recover_x) recover()

    if(tune) mixture_weight_trace[[1]][[1]] <<- c(mixture_weight_trace[[1]][[1]][-1],0)

    if(log(runif(1)) < logR){

      weight_active <<- weight_active_proposed
      weight_within_active <<- weight_within_active_proposed

      if(tune) mixture_weight_trace[[1]][[1]][100] <<- 1

    }

  },

  mhCombinedSigmaMu = function(recover_x = FALSE){

    lib = 1

    active_means_dif_proposed <<- active_means_dif

    inactive_means_proposed <<- inactive_means

    inactive_variances_proposed <<- inactive_variances

    active_variances_proposed <<- active_variances

    not_inspike <- which(inactive_spike_allocation == 0)

    random_mu <- sample(0:num_active_components, 1)

    random_sigma <- sample(0:num_active_components, 1)

    HR = 0


    ### propose new mean

    if(random_mu == 0){

      inactive_means_proposed <<- threshold_i - abs(threshold_i - inactive_means - rnorm(1, 0, inactive_mean_tuningParam))

      mu_pp = computeInactiveMeansPriorProbability(inactive_means_proposed)

      mu_pc = computeInactiveMeansPriorProbability(inactive_means)

    }else{

      active_means_dif_proposed[random_mu] <<- abs(active_means_dif[random_mu] + rnorm(1, 0, active_mean_tuningParam[random_mu]))

      mu_pp = computeActiveMeansDifPriorProbability(active_means_dif_proposed)

      mu_pc = computeActiveMeansDifPriorProbability(active_means_dif)

    }

    ## propose new sigma

    if(random_sigma == 0){

      inactive_variances_proposed <<- 10^(.self$two_boundary_slide_move(log(inactive_variances, 10), inactive_variances_prior_log_min, inactive_variances_prior_log_max, inactive_variance_tuningParam))

      sigma_pp = computeInactiveVariancesPriorProbability(inactive_variances_proposed)

      sigma_pc = computeInactiveVariancesPriorProbability(inactive_variances)

      HR = log(inactive_variances_proposed/inactive_variances)


    }else{


      active_variances_proposed[random_sigma] <<- 10^(.self$two_boundary_slide_move(log(active_variances[random_sigma], 10), active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[random_sigma]))

      sigma_pp = computeActiveVariancesPriorProbability(active_variances_proposed)

      sigma_pc = computeActiveVariancesPriorProbability(active_variances)

      HR = log(active_variances_proposed[random_sigma]/active_variances[random_sigma])


    }


    ## inactive component likelihood

    Lp = (1 - weight_active) * (1 - spike_probability) *
      1/(sqrt(inactive_variances_proposed) * sqrt2pi * exp((Yg[not_inspike] - inactive_means_proposed)^2/(2 * inactive_variances_proposed)))

    Lc = (1 - weight_active) * (1 - spike_probability) *
      1/(sqrt(inactive_variances) * sqrt2pi * exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances)))


    ## active component likelihood

    mp = .self$calculate_active_means(active_means_dif_proposed)
    mc = active_means

    Lpk = 0
    Lck = 0

    if(length(active_idx) > 0){

      for(k in 1:num_active_components){

        Lpk = Lpk + weight_within_active[k] *
          1/(sqrt(active_variances_proposed[k]) * sqrt2pi * exp((Yg[not_inspike] - mp[k])^2/(2 * active_variances_proposed[k])))

        Lck = Lck + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - mc[k])^2/(2 * active_variances[k])))

      }

    }

    lnLp <- beta_oneLib * sum(log(Lp + weight_active * Lpk))
    lnLc <- beta_oneLib * sum(log(Lc + weight_active * Lck))

    R = temperature * (lnLp + mu_pp + sigma_pp - lnLc - mu_pc - sigma_pc) + HR

    if(recover_x) recover()

    if(log(runif(1)) < R){

      # cat("*********SIIGMAA_____MUUUUUUU******", random_mu, "  ", random_sigma, "\n")


      if(random_mu == 0 && random_sigma == 0){


        inactive_means <<- inactive_means_proposed


        inactive_variances <<- inactive_variances_proposed


      }else if(random_mu == 0 && random_sigma > 0){


        inactive_means <<- inactive_means_proposed


        active_variances[random_sigma] <<- active_variances_proposed[random_sigma]



      }else if( random_mu > 0 && random_sigma == 0){


        active_means_dif[random_mu] <<- active_means_dif_proposed[random_mu]

        active_means <<- mp


        inactive_variances <<- inactive_variances_proposed


      }else{


        active_means_dif[random_mu] <<- active_means_dif_proposed[random_mu]

        active_means <<- mp


        active_variances[random_sigma] <<- active_variances_proposed[random_sigma]

      }

      .self$gibbsAllocationActiveInactive()
      .self$gibbsAllocationWithinActive()

    }

  },

  mhSpikeAllocation_oneLibrary = function(recover_x = FALSE, tune = FALSE){

    # note, proposed moves are only within the inactive allocation,
    # no moves between the spike and the active component are ever proposed

    zero_and_inactive_idx = all_zero_idx[which(allocation_active_inactive[all_zero_idx] == 0)]

    inactive_spike_allocation_proposal = abs(inactive_spike_allocation[zero_and_inactive_idx] - 1)

    proposed_out_spike_idx = which(inactive_spike_allocation_proposal==0)
    proposed_in_spike_idx = which(inactive_spike_allocation_proposal==1)

    current_out_spike_idx = which(inactive_spike_allocation[zero_and_inactive_idx]==0)
    current_in_spike_idx = which(inactive_spike_allocation[zero_and_inactive_idx]==1)


    #*** Make proposals

    yg_proposed = rnorm(length(zero_and_inactive_idx), inactive_means, sqrt(inactive_variances))

    p_x_proposed <- .self$get_px(yy = yg_proposed, gl = gene_lengths[zero_and_inactive_idx])
    # p_x_proposed <- .self$xxget_px(yy = yg_proposed, gl = gene_lengths[zero_and_inactive_idx], xx = Xg[zero_and_inactive_idx,], sg = sigma_g[zero_and_inactive_idx])


    #*** Compute likelihoods and prior probabilities

    p_yg_proposed <- .self$computeInactiveLikelihood(yg_proposed, spike_probability, inactive_means, inactive_variances,
                                                     inactive_spike_allocation_proposal)

    p_yg <- .self$computeInactiveLikelihood(Yg[zero_and_inactive_idx], spike_probability, inactive_means, inactive_variances,
                                            inactive_spike_allocation[zero_and_inactive_idx])


    pp = NULL; pc = NULL

    pp[c(proposed_out_spike_idx, proposed_in_spike_idx)] = c(p_yg_proposed[proposed_out_spike_idx] ,
                                                             rep(log(spike_probability), length(proposed_in_spike_idx)))

    pc[c(current_in_spike_idx, current_out_spike_idx)] = c(rep(log(spike_probability), length(current_in_spike_idx)),
                                                           p_yg[current_out_spike_idx] )


    #/// Compute Hastings Ratio

    HR = NULL

    HR[c(proposed_out_spike_idx, proposed_in_spike_idx)] = c(-(p_yg_proposed[proposed_out_spike_idx]  - log(1-spike_probability)),
                                                             p_yg[proposed_in_spike_idx]  - log(1-spike_probability))



    Lp <- .self$computeXgLikelihood_oneLibrary(Xg[zero_and_inactive_idx], p_x_proposed, spike_allocation = inactive_spike_allocation_proposal,
                                               nd_spike = no_detect_spike[zero_and_inactive_idx])

    Lc <- .self$computeXgLikelihood_oneLibrary(Xg[zero_and_inactive_idx], p_x[zero_and_inactive_idx], spike_allocation = inactive_spike_allocation[zero_and_inactive_idx],
                                               nd_spike = no_detect_spike[zero_and_inactive_idx])



    #*** Update Parameters

    R = Lp + pp - Lc - pc + HR

    choice <- (log(runif(length(zero_and_inactive_idx))) < R) + 1

    choice_matrix <- cbind(seq(length(zero_and_inactive_idx)), choice)

    current_proposed = cbind(inactive_spike_allocation[zero_and_inactive_idx], inactive_spike_allocation_proposal)

    yg_current_proposed = cbind(Yg[zero_and_inactive_idx], yg_proposed)


    if(recover_x) recover()


    inactive_spike_allocation[zero_and_inactive_idx] <<- current_proposed[choice_matrix]

    Yg[zero_and_inactive_idx][proposed_out_spike_idx] <<- yg_current_proposed[choice_matrix][proposed_out_spike_idx]

    Yg[zero_and_inactive_idx][proposed_in_spike_idx] <<- yg_current_proposed[choice_matrix][proposed_in_spike_idx]

    .self$set_sigmaX_pX()

    .self$setInSpike_idx()


  },

  mhP_x_oneLibrary = function(recover_x = FALSE, tune = FALSE){

    randlib=sample(num_libraries,1)

    proposed_alpha_r <- alpha_r

    proposal <- .self$scale_move(alpha_r[randlib], tuningParam_alpha_r[randlib])

    proposed_alpha_r[randlib] <- proposal[[1]]

    if(tune) alpha_r_trace[[1]][[1]][randlib,] <<- c(alpha_r_trace[[1]][[1]][randlib, -1], 0)

    HR <- log(proposal[[2]])

    proposed_p_x <- .self$get_px(falpha_r = proposed_alpha_r, yy = Yg)
    # proposed_p_x <- .self$xxget_px(falpha_r = proposed_alpha_r, yy = Yg, xx = Xg[,randlib])

    Lp <- sum(.self$computeXgLikelihood_oneLibrary(Xg, proposed_p_x))

    Lc <- sum(.self$computeXgLikelihood_oneLibrary(Xg, p_x))

    pp <- .self$computeAlphaRPriorProbability(proposed_alpha_r)

    pc <- .self$computeAlphaRPriorProbability(alpha_r)

    R <- Lp - Lc + pp - pc + HR

    if(recover_x) recover()

    if(log(runif(1)) < R & ! is.nan(R) & ! is.na(R)){

      alpha_r[randlib] <<- proposed_alpha_r[randlib]

      p_x[,randlib] <<- proposed_p_x[,randlib]

      if(tune) alpha_r_trace[[1]][[1]][randlib,100] <<- 1

    }

  },

  mhYg_oneLibrary = function(recover_x = FALSE, tune = FALSE){

    nodetect_inactive_idx <- which(allocation_active_inactive[all_zero_idx] == 0)

    nodetect_active_idx <- which(allocation_active_inactive[all_zero_idx] == 1)

    Yg_proposed <<- Yg[all_zero_idx] + rnorm(length(all_zero_idx), 0, tuningParam_yg[all_zero_idx]) * (1 - inactive_spike_allocation[all_zero_idx])

    p_x_proposed <- .self$get_px(yy = Yg_proposed,  gl = gene_lengths[all_zero_idx])
    # p_x_proposed <- .self$xxget_px(yy = Yg_proposed,  gl = gene_lengths[all_zero_idx], xx = Xg[all_zero_idx,], sg = sigma_g[all_zero_idx])

    ## proposed Yg likelihoods
    proposed_XgLikelihood = .self$computeXgLikelihood_oneLibrary(Xg[all_zero_idx], p_x_proposed, spike_allocation = inactive_spike_allocation[all_zero_idx],
                                                                 spike_prob = spike_probability, nd_spike = no_detect_spike[all_zero_idx])

    Lp_pp = proposed_XgLikelihood

    Lp_pp[nodetect_inactive_idx] = Lp_pp[nodetect_inactive_idx] +
      .self$computeInactiveLikelihood(Yg_proposed[nodetect_inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[all_zero_idx][nodetect_inactive_idx])

    Lp_pp[nodetect_active_idx] = Lp_pp[nodetect_active_idx] +
      .self$computeActiveLikelihood(Yg_proposed[nodetect_active_idx],  allocation_within_active[[1]][all_zero_idx][nodetect_active_idx], active_means, active_variances)


    ## current Yg likelihoods
    current_XgLikelihood <- .self$computeXgLikelihood_oneLibrary(Xg[all_zero_idx], p_x[all_zero_idx], spike_allocation = inactive_spike_allocation[all_zero_idx],
                                                                 spike_prob = spike_probability, nd_spike = no_detect_spike[all_zero_idx])

    Lc_pc <- current_XgLikelihood

    Lc_pc[nodetect_inactive_idx] = Lc_pc[nodetect_inactive_idx] +
      .self$computeInactiveLikelihood(Yg[all_zero_idx][nodetect_inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[all_zero_idx][nodetect_inactive_idx])

    Lc_pc[nodetect_active_idx] = Lc_pc[nodetect_active_idx] +
      .self$computeActiveLikelihood(Yg[all_zero_idx][nodetect_active_idx],  allocation_within_active[[1]][all_zero_idx][nodetect_active_idx], active_means, active_variances)

    Lp_pp[which(is.na(Lp_pp))] <- -Inf

    R = Lp_pp - Lc_pc

    current_and_proposed <- cbind(Yg[all_zero_idx], Yg_proposed)

    current_and_proposed_likelihood <- cbind(current_XgLikelihood, proposed_XgLikelihood)

    if(recover_x) recover()

    choice_idx <- 1 + (log(runif(length(Yg_proposed))) < R)

    Yg[all_zero_idx] <<- current_and_proposed[cbind(seq(length(Yg_proposed)), choice_idx)]


    ## update XgLikelihood, Sg and p_x

    XgLikelihood[all_zero_idx] <<- current_and_proposed_likelihood[cbind(seq(length(Yg_proposed)), choice_idx)]

    .self$set_sigmaX_pX()

    if(tune) Yg_trace[all_zero_idx,] <<- cbind(Yg_trace[all_zero_idx,-1], choice_idx - 1)

  },

  mhCombinedMultiSigma = function(sigma_test = FALSE, tune = FALSE){

    inactive_variances_proposed <<- inactive_variances

    active_variances_proposed <<- active_variances

    not_inspike <- which(inactive_spike_allocation == 0)

    HR = 0

    random_components = sample(0:num_active_components, 2, replace = FALSE)

    if(0 %in% random_components){

      pk <- random_components[which(random_components != 0)]

      if(rbinom(1, 1, 0.5) < 1){

        proposal <- .self$twoSigma_slide_move(log(c(inactive_variances, active_variances_proposed[pk]), 10), inactive_variances_prior_log_min, inactive_variances_prior_log_max, inactive_variance_tuningParam)

        inactive_variances_proposed <<- 10^proposal[1]

        active_variances_proposed[pk] <<- 10^proposal[2]

      }else{

        proposal <- .self$twoSigma_slide_move(log(c(active_variances_proposed[pk], inactive_variances), 10), active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[pk])

        inactive_variances_proposed <<- 10^proposal[2]

        active_variances_proposed[pk] <<- 10^proposal[1]

      }

      pp = computeInactiveVariancesPriorProbability(inactive_variances_proposed) + computeActiveVariancesPriorProbability(active_variances_proposed)

      pc = computeInactiveVariancesPriorProbability(inactive_variances) + computeActiveVariancesPriorProbability(active_variances)

    }else{

      pk1 <- random_components[1]

      pk2 <- random_components[2]

      proposal <- .self$twoSigma_slide_move(log(active_variances_proposed[c(pk1, pk2)], 10), active_variances_prior_log_min, active_variances_prior_log_max, active_variance_tuningParam[pk1])

      active_variances_proposed[pk1] <<- 10^proposal[1]

      active_variances_proposed[pk2] <<- 10^proposal[2]

      pp = computeActiveVariancesPriorProbability(active_variances_proposed)

      pc = computeActiveVariancesPriorProbability(active_variances)

    }


    HR = sum(log(c(inactive_variances_proposed, active_variances_proposed)/c(inactive_variances, active_variances)))

    if(tune) multi_sigma_trace[[1]][[1]] <<- c(multi_sigma_trace[[1]][[1]][-1],0)


    Lp = (1 - spike_probability) * 1/(sqrt(inactive_variances_proposed) * sqrt2pi *
                                        exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances_proposed)))

    Lc = (1 - spike_probability) * 1/(sqrt(inactive_variances) * sqrt2pi *
                                        exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances)))

    if(!is.null(active_gene_set)){

      actg_idx = which(gene_names[out_spike_idx] %in% active_gene_set)

      Lp[actg_idx] = 0

      Lc[actg_idx] = 0

    }

    Lpk = 0
    Lck = 0

    if(length(active_idx) > 0){

      for(k in 1:num_active_components){

        Lpk = Lpk + weight_within_active[k] *
          1/(sqrt(active_variances_proposed[k]) * sqrt2pi * exp((Yg[not_inspike] - active_means[k])^2/(2 * active_variances_proposed[k])))

        Lck = Lck + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - active_means[k])^2/(2 * active_variances[k])))

      }

    }

    lnLp <- beta_oneLib * sum(log((1 - weight_active) * Lp + weight_active * Lpk))
    lnLc <- beta_oneLib * sum(log((1 - weight_active) * Lc + weight_active * Lck))

    if(sigma_test) recover()

    R = temperature * (lnLp + pp - lnLc - pc) + HR


    if(log(runif(1)) < R & ! is.nan(R) & ! is.na(R)){

      # cat("*************MULTISIIIIIIGMAAAAAA!*********", random_components, "\n")

      if(0 %in% random_components){

        inactive_variances <<- inactive_variances_proposed

        active_variances[pk] <<- active_variances_proposed[pk]

      }else{

        active_variances[pk1] <<- active_variances_proposed[pk1]

        active_variances[pk2] <<- active_variances_proposed[pk2]

      }

      if(tune) multi_sigma_trace[[1]][[1]][100] <<- 1

      .self$gibbsAllocationActiveInactive()
      .self$gibbsAllocationWithinActive()


    }

  },

  mhCombinedMultiMean = function(mean_test = FALSE, tune = FALSE){

    # propose new component parameters
    # propose/accept new allocations given the new parameters
    # calculate R to accept/reject the new component params and new allocations

    lib <- 1

    active_means_dif_proposed <<- active_means_dif

    inactive_means_proposed <<- inactive_means

    not_inspike <- which(inactive_spike_allocation == 0)

    # inspike <- which(inactive_spike_all == 1)


    pk <- sample(seq(num_active_components), 1)

    if(pk < num_active_components){

      active_means_dif_proposed[c(pk, pk + 1)] <<- .self$twoMean_slide_move(active_means_dif[c(pk, pk + 1)], active_mean_tuningParam[pk])

      if(tune) active_means_trace[[1]][[lib]][pk,] <<- c(active_means_trace[[1]][[lib]][pk,-1],0)

    }else{

      active_means_dif_proposed[pk] <<- abs(active_means_dif[pk] + rnorm(1,0,active_mean_tuningParam[pk]))

      if(tune) active_means_trace[[1]][[lib]][pk,] <<- c(active_means_trace[[1]][[lib]][pk,-1],0)

    }

    pp = computeActiveMeansDifPriorProbability(active_means_dif_proposed)
    pc = computeActiveMeansDifPriorProbability(active_means_dif)

    Lp = (1 - weight_active) * (1 - spike_probability) *
      1/(sqrt(inactive_variances) * sqrt2pi * exp((Yg[not_inspike] - inactive_means_proposed)^2/(2 * inactive_variances)))

    Lc = (1 - weight_active) * (1 - spike_probability) *
      1/(sqrt(inactive_variances) * sqrt2pi * exp((Yg[not_inspike] - inactive_means)^2/(2 * inactive_variances)))


    if(!is.null(active_gene_set)){

      actg_idx = which(gene_names[out_spike_idx] %in% active_gene_set)

      Lp[actg_idx] = 0

      Lc[actg_idx] = 0

    }

    Lpk = 0
    Lck = 0

    if(length(active_idx) > 0){

      mp = .self$calculate_active_means(active_means_dif_proposed)
      mc = active_means

      for(k in 1:num_active_components){

        Lpk = Lpk + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - mp[k])^2/(2 * active_variances[k])))

        Lck = Lck + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[not_inspike] - mc[k])^2/(2 * active_variances[k])))

      }

    }

    lnLp <- beta_oneLib * sum(log(Lp + weight_active * Lpk))
    lnLc <- beta_oneLib * sum(log(Lc + weight_active * Lck))

    R = temperature * (lnLp + pp - lnLc - pc)

    if(mean_test) recover()


    if(log(runif(1)) < R & ! is.nan(R) & ! is.na(R)){

      # print("#______________________________ OLD _________________________________#")

      active_means_dif <<- active_means_dif_proposed

      if(tune) active_means_trace[[1]][[lib]][pk,100] <<- 1

      if(length(active_idx) > 0) active_means <<- mp

      .self$gibbsAllocationActiveInactive()
      .self$gibbsAllocationWithinActive()


    }


  },

  mhSigmaSwap = function(ss_test = FALSE, tune = FALSE){

    inactive_variances_proposed <<- inactive_variances

    active_variances_proposed <<- active_variances

    weight_within_active_proposed <<- weight_within_active

    not_inspike <- which(inactive_spike_allocation == 0)

    HR = 0

    random_components = sample(1:num_active_components, 2, replace = FALSE)

    pk1 <- random_components[1]

    pk2 <- random_components[2]

    weight_within_active_proposed[pk1] <<- weight_within_active[pk2]

    weight_within_active_proposed[pk2] <<- weight_within_active[pk1]

    active_variances_proposed[pk1] <<- active_variances[pk2]

    active_variances_proposed[pk2] <<- active_variances[pk1]

    pp = computeActiveVariancesPriorProbability(active_variances_proposed)

    pc = computeActiveVariancesPriorProbability(active_variances)


    Lpk = 0
    Lck = 0

    if(length(active_idx) > 0){

      for(k in 1:num_active_components){

        Lpk = Lpk + weight_within_active_proposed[k] *
          1/(sqrt(active_variances_proposed[k]) * sqrt2pi * exp((Yg[active_idx] - active_means[k])^2/(2 * active_variances_proposed[k])))

        Lck = Lck + weight_within_active[k] *
          1/(sqrt(active_variances[k]) * sqrt2pi * exp((Yg[active_idx] - active_means[k])^2/(2 * active_variances[k])))

      }

    }

    lnLp <- beta_oneLib * sum(log(Lpk))
    lnLc <- beta_oneLib * sum(log(Lck))

    R = temperature * (lnLp + pp - lnLc - pc) + HR

    if(ss_test) recover()

    if(log(runif(1)) < R & ! is.nan(R) & ! is.na(R)){

      cat("************SWAAAAAPPPP!!!********", active_means[random_components], "\n")

      weight_within_active <<- weight_within_active_proposed

      active_variances[pk1] <<- active_variances_proposed[pk1]

      active_variances[pk2] <<- active_variances_proposed[pk2]

      .self$gibbsAllocationWithinActive()

    }


  }

)













