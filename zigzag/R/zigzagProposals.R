zigzag$methods(

  ##########################
  ### Lower level params ###
  ##########################

  mhVariance_g = function(recover_x = FALSE, tune = FALSE){

    HR = 0

    proposal <- .self$scale_move(variance_g[out_spike_idx], tuningParam_variance_g[out_spike_idx])

    variance_g_proposed <- proposal[[1]]

    HR <- log(proposal[[2]])

    proposed_XgLikelihood <- .self$computeXgLikelihood(Xg[out_spike_idx,], Yg[out_spike_idx], variance_g_proposed, p_x[out_spike_idx,],
                                                       spike_allocation = inactive_spike_allocation[out_spike_idx], nd_spike = no_detect_spike[out_spike_idx])

    proposed_varianceGProbability <- .self$computeVarianceGPriorProbability(variance_g_proposed, Sg[out_spike_idx])

    Lp_pp = proposed_XgLikelihood + proposed_varianceGProbability

    Lc_pp = XgLikelihood[out_spike_idx] + variance_g_probability[out_spike_idx]

    R = temperature * (Lp_pp - Lc_pp) + HR

    R[which(variance_g_proposed > variance_g_upper_bound)] <- -Inf

    if(recover_x) recover()

    current_and_proposed <- cbind(variance_g[out_spike_idx], variance_g_proposed)

    current_and_proposed_likelihood <- cbind(XgLikelihood[out_spike_idx], proposed_XgLikelihood)


    num_choices <- length(out_spike_idx)

    choice_idx <- (log(runif(num_choices)) < R) + 1  # 1 means reject proposal, 2 means accept proposal.

    choice_matrix <- cbind(seq(num_choices), choice_idx) # row-column pairs for current_and_proposed_variance_g_probability with 1st column meaning reject proposal

    variance_g[out_spike_idx] <<- current_and_proposed[choice_matrix]

    if(tune){

      if(length(out_spike_idx) < 2){

        variance_g_trace[out_spike_idx,] <<- cbind(matrix(variance_g_trace[out_spike_idx,-1], nrow = length(out_spike_idx)), choice_idx - 1)

      }else{

        variance_g_trace[out_spike_idx,] <<- cbind(variance_g_trace[out_spike_idx,-1], choice_idx - 1)
      }
    }


    XgLikelihood[out_spike_idx] <<- current_and_proposed_likelihood[choice_matrix]

    current_and_proposed_variance_g_probability <- cbind(variance_g_probability[out_spike_idx], proposed_varianceGProbability)

    variance_g_probability[out_spike_idx] <<- current_and_proposed_variance_g_probability[choice_matrix]

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

    proposed_variance_g_Likelihood <- .self$computeVarianceGPriorProbability(variance_g[out_spike_idx], proposed_Sg[out_spike_idx])

    Lp <- sum(proposed_variance_g_Likelihood)

    Lc <- sum(variance_g_probability[out_spike_idx])

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

      variance_g_probability[out_spike_idx] <<- proposed_variance_g_Likelihood

    }

  },

  mhTau = function(recover_x = FALSE, tune = FALSE){

    proposal <- .self$scale_move(tau, tuningParam_tau)
    proposed_tau <- proposal[[1]]

    HR <- log(proposal[[2]])

    proposed_tau_Likelihood <- .self$computeVarianceGPriorProbability(variance_g[out_spike_idx], Sg[out_spike_idx], gtau = proposed_tau)
    tau_Likelihood <- variance_g_probability[out_spike_idx]

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

      variance_g_probability[out_spike_idx] <<- proposed_tau_Likelihood

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

    proposed_variance_g_Likelihood <- .self$computeVarianceGPriorProbability(variance_g[out_spike_idx], proposed_Sg[out_spike_idx], proposed_tau)

    Lp <- sum(proposed_variance_g_Likelihood)

    Lc <- sum(variance_g_probability[out_spike_idx])

    pp <- .self$computeS0PriorProbability(proposed_s0) + .self$computeTauPriorProbability(proposed_tau)

    pc <- .self$computeS0PriorProbability(s0) + .self$computeTauPriorProbability(tau)

    R <- temperature * (Lp - Lc + pp - pc) + HR

    if(recover_x) recover()

    if(log(runif(1)) < R){

      s0 <<- proposed_s0

      tau <<- proposed_tau

      Sg <<- proposed_Sg

      if(tune) s0tau_trace[[1]][[1]][100] <<- 1

      variance_g_probability[out_spike_idx] <<- proposed_variance_g_Likelihood

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

      proposed_xgLikelihood <- XgLikelihood + .self$computeXgLikelihood(Xg[,randlib], Yg, variance_g, proposed_p_x[,randlib]) -
        .self$computeXgLikelihood(Xg[,randlib], Yg, variance_g, p_x[,randlib])

    }else{

      proposed_xgLikelihood <- .self$computeXgLikelihood(Xg, Yg, variance_g, proposed_p_x)

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

    num_choice = length(out_spike_idx)

    inactive_outspike_idx = which(allocation_active_inactive[out_spike_idx] == 0)
    active_outspike_idx = which(allocation_active_inactive[out_spike_idx] == 1)

    Yg_proposed <<-  rnorm(num_choice, Yg[out_spike_idx], tuningParam_yg[out_spike_idx])#only change out of spike genes

    Sg_proposed <- .self$get_sigmax(yy = Yg_proposed)

    p_x_proposed <- .self$get_px(yy = Yg_proposed, gl = gene_lengths[out_spike_idx])

    ## proposed Yg likelihoods
    proposed_XgLikelihood = .self$computeXgLikelihood(Xg[out_spike_idx,], Yg_proposed, variance_g[out_spike_idx], p_x_proposed,
                                                      spike_allocation = inactive_spike_allocation[out_spike_idx])

    proposed_varianceGProbability = .self$computeVarianceGPriorProbability(variance_g[out_spike_idx], Sg_proposed)

    Lp_pp = proposed_XgLikelihood

    Lp_pp[inactive_outspike_idx] = Lp_pp[inactive_outspike_idx] +
      .self$computeInactiveLikelihood(Yg_proposed[inactive_outspike_idx], spike_probability, inactive_means,
                                      inactive_variances, inactive_spike_allocation[out_spike_idx][inactive_outspike_idx])

    Lp_pp[active_outspike_idx] = Lp_pp[active_outspike_idx] +
      .self$computeActiveLikelihood(Yg_proposed[active_outspike_idx],  allocation_within_active[[1]][out_spike_idx][active_outspike_idx],
                                    active_means, active_variances)

    Lp_pp = Lp_pp + proposed_varianceGProbability


    ## current Yg likelihoods
    Lc_pc = XgLikelihood[out_spike_idx]

    Lc_pc[inactive_outspike_idx] = Lc_pc[inactive_outspike_idx] +
      .self$computeInactiveLikelihood(Yg[out_spike_idx][inactive_outspike_idx], spike_probability, inactive_means,
                                      inactive_variances, inactive_spike_allocation[out_spike_idx][inactive_outspike_idx])

    Lc_pc[active_outspike_idx] = Lc_pc[active_outspike_idx] +
      .self$computeActiveLikelihood(Yg[out_spike_idx][active_outspike_idx],  allocation_within_active[[1]][out_spike_idx][active_outspike_idx],
                                    active_means, active_variances)

    Lc_pc = Lc_pc + variance_g_probability[out_spike_idx]

    R = temperature * (Lp_pp - Lc_pc)

    current_and_proposed <- cbind(Yg[out_spike_idx], Yg_proposed)

    current_and_proposed_likelihood <- cbind(XgLikelihood[out_spike_idx], proposed_XgLikelihood)

    if(recover_x) recover()


    choice_idx <- 1 + (log(runif(num_choice)) < R)

    choice_matrix <- cbind(seq(num_choice), choice_idx)

    Yg[out_spike_idx] <<- current_and_proposed[cbind(seq(num_choice), choice_idx)]

    if(tune){

      if(length(out_spike_idx) < 2){

        Yg_trace[out_spike_idx,] <<- cbind(matrix(Yg_trace[out_spike_idx,-1], nrow = length(out_spike_idx)), choice_idx - 1)

      }else{

        Yg_trace[out_spike_idx,] <<- cbind(Yg_trace[out_spike_idx,-1], choice_idx - 1)
      }
    }

    ## update XgLikelihood, Sg and p_x

    XgLikelihood[out_spike_idx] <<- current_and_proposed_likelihood[choice_matrix]

    .self$set_sigmaX_pX()

    current_and_proposed_varianceGProbability <- cbind(variance_g_probability[out_spike_idx], proposed_varianceGProbability)

    variance_g_probability[out_spike_idx] <<- current_and_proposed_varianceGProbability[choice_matrix]

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

    variance_g_proposed <- rlnorm(length(zero_and_inactive_idx), log(Sg_proposed) + tau, sqrt(tau))  # P(variance_g | proposed_Yg)

    p_x_proposed <- .self$get_px(yy = yg_proposed, gl = gene_lengths[zero_and_inactive_idx])


    #*** Compute prior probabilities

    p_variance_g_proposed <- .self$computeVarianceGPriorProbability(variance_g_proposed, Sg_proposed)

    p_variance_g <- .self$computeVarianceGPriorProbability(variance_g[zero_and_inactive_idx], Sg[zero_and_inactive_idx])

    p_yg_proposed <- .self$computeInactiveLikelihood(yg_proposed, spike_probability, inactive_means,
                                                     inactive_variances, inactive_spike_allocation_proposal)

    p_yg <- .self$computeInactiveLikelihood(Yg[zero_and_inactive_idx], spike_probability, inactive_means,
                                            inactive_variances, inactive_spike_allocation[zero_and_inactive_idx])

    pp = NULL; pc = NULL

    ### this needs a clear explanation.
    pp[c(proposed_out_spike_idx, proposed_in_spike_idx)] = c(p_yg_proposed[proposed_out_spike_idx] + p_variance_g_proposed[proposed_out_spike_idx],
                                                             rep(log(spike_probability), length(proposed_in_spike_idx)))

    pc[c(current_in_spike_idx, current_out_spike_idx)] = c(rep(log(spike_probability), length(current_in_spike_idx)),
                                                           p_yg[current_out_spike_idx] + p_variance_g[current_out_spike_idx])



    #*** Compute Likelihoods

    Lp = .self$computeXgLikelihood(Xg[zero_and_inactive_idx,], yg_proposed, variance_g_proposed, p_x_proposed,
                                   inactive_spike_allocation_proposal, nd_spike = no_detect_spike[zero_and_inactive_idx])
    Lc = XgLikelihood[zero_and_inactive_idx]

    both_inf_boolean = ifelse(Lp == -Inf & Lc == -Inf, TRUE, FALSE)
    Lp[both_inf_boolean] = 0
    Lc[both_inf_boolean] = 0

    #/// Compute Hastings Ratio

    HR = NULL

    # subtract out spike_probability because of use in computeInactiveLikelihood
    HR[c(proposed_out_spike_idx, proposed_in_spike_idx)] = c(-(p_yg_proposed[proposed_out_spike_idx] - log(1-spike_probability) +
                                                                 p_variance_g_proposed[proposed_out_spike_idx]),
                                                             p_yg[proposed_in_spike_idx] - log(1-spike_probability) +
                                                               p_variance_g[proposed_in_spike_idx])



    #*** Update Parameters

    R <- temperature * (Lp - Lc + pp - pc) + HR

    R[which(variance_g_proposed > variance_g_upper_bound)] <- -Inf

    choice <- (log(runif(length(zero_and_inactive_idx))) < R) + 1

    # choice[which(is.na(choice) | is.nan(choice))] <- 1 #Do not update (set choice to current value) if NA or NaN.

    choice_matrix <- cbind(seq(length(zero_and_inactive_idx)), choice)

    current_proposed = cbind(inactive_spike_allocation[zero_and_inactive_idx], inactive_spike_allocation_proposal)

    if(recover_x) recover()


    inactive_spike_allocation[zero_and_inactive_idx] <<- current_proposed[choice_matrix]

    yg_current_proposed = cbind(Yg[zero_and_inactive_idx], yg_proposed)

    Yg[zero_and_inactive_idx][proposed_out_spike_idx] <<- yg_current_proposed[choice_matrix][proposed_out_spike_idx]

    variance_g_current_proposed <- cbind(variance_g[zero_and_inactive_idx], variance_g_proposed)

    variance_g[zero_and_inactive_idx][proposed_out_spike_idx] <<- variance_g_current_proposed[choice_matrix][proposed_out_spike_idx]

    .self$set_sigmaX_pX()

    .self$setInSpike_idx()

    xgl_current_proposed = cbind(Lc, Lp)

    XgLikelihood[zero_and_inactive_idx] <<- xgl_current_proposed[choice_matrix]

    variance_gProb_current_proposed <- cbind(p_variance_g, p_variance_g_proposed)

    variance_g_probability[zero_and_inactive_idx][proposed_out_spike_idx] <<- variance_gProb_current_proposed[choice_matrix][proposed_out_spike_idx]

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

    k = sample(num_active_components, 1, replace = TRUE)

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


  ###################
  ### DEVELOPMENT ###
  ###################

  mhS0Tau_shift = function(recover_x = FALSE, tune = FALSE){
    HR = 0

    # proposal <- .self$scale_move(tau, tuningParam_s0tau)
    #
    # proposed_s0 <- s0 - log(proposal[[2]])
    #
    # proposed_tau <- proposal[[1]]
    #
    # HR <- log(proposal[[2]])

    shift <- rnorm(1, 0, tuningParam_s0tau)

    proposed_s0 <- s0 + shift

    proposed_tau <- abs(tau - 0.63 * shift)

    if(tune) s0tau_trace[[1]][[1]] <<- c(s0tau_trace[[1]][[1]][-1], 0)

    proposed_Sg <- .self$get_sigmax(ss0 = proposed_s0, ss1 = s1)

    proposed_variance_g_Likelihood <- .self$computeVarianceGPriorProbability(variance_g[out_spike_idx], proposed_Sg[out_spike_idx], proposed_tau)

    Lp <- sum(proposed_variance_g_Likelihood)

    Lc <- sum(variance_g_probability[out_spike_idx])

    pp <- .self$computeS0PriorProbability(proposed_s0) + .self$computeTauPriorProbability(proposed_tau)

    pc <- .self$computeS0PriorProbability(s0) + .self$computeTauPriorProbability(tau)

    R <- temperature * (Lp - Lc + pp - pc) + HR

    if(recover_x) recover()

    if(log(runif(1)) < R){

      s0 <<- proposed_s0

      tau <<- proposed_tau

      Sg <<- proposed_Sg

      if(tune) s0tau_trace[[1]][[1]][100] <<- 1

      variance_g_probability[out_spike_idx] <<- proposed_variance_g_Likelihood

    }

  },

  mhYgVarianceg = function(recover_x = FALSE, tune = FALSE){

    # mean_Xg = rowMeans(Xg[out_spike_idx,])
    # sd_Xg = sapply(out_spike_idx, function(g){sd(Xg[g,])}) #alternatively, compute sd_meanXg

    Sg_current = .self$get_sigmax(yy = Yg[out_spike_idx])

    # proposed new yg
    num_choice = length(out_spike_idx)

    inactive_outspike_idx = which(allocation_active_inactive[out_spike_idx] == 0)
    active_outspike_idx = which(allocation_active_inactive[out_spike_idx] == 1)

    Yg_proposed <<-  rnorm(num_choice, rwm[out_spike_idx], sqrt(rwv[out_spike_idx]))#only change out of spike genes

    Sg_proposed <- .self$get_sigmax(yy = Yg_proposed)

    p_x_proposed <- .self$get_px(yy = Yg_proposed, gl = gene_lengths[out_spike_idx])


    # propose new variance_g
    variance_g_proposed <- rlnorm(num_choice, log(Sg_proposed) + tau, sqrt(tau))


    # compute hastings ratio
    HR <- dlnorm(variance_g[out_spike_idx], log(Sg_current) + tau, sqrt(tau), log = TRUE) +
      dnorm(Yg[out_spike_idx], rwm[out_spike_idx], sqrt(rwv[out_spike_idx]), log = TRUE) -
      dlnorm(variance_g_proposed, log(Sg_proposed) + tau, sqrt(tau), log = TRUE) -
      dnorm(Yg_proposed, rwm[out_spike_idx], sqrt(rwv[out_spike_idx]), log = TRUE)


    # compute current priors X Likelihoods
    Lc_pc = XgLikelihood[out_spike_idx]

    Lc_pc[inactive_outspike_idx] = Lc_pc[inactive_outspike_idx] +
      .self$computeInactiveLikelihood(Yg[out_spike_idx][inactive_outspike_idx], spike_probability, inactive_means,
                                      inactive_variances, inactive_spike_allocation[out_spike_idx][inactive_outspike_idx])

    Lc_pc[active_outspike_idx] = Lc_pc[active_outspike_idx] +
      .self$computeActiveLikelihood(Yg[out_spike_idx][active_outspike_idx],  allocation_within_active[[1]][out_spike_idx][active_outspike_idx],
                                    active_means, active_variances)

    Lc_pc = Lc_pc + variance_g_probability[out_spike_idx]



    # compute proposed priors X Likelihoods
    proposed_XgLikelihood <- .self$computeXgLikelihood(Xg[out_spike_idx,], Yg_proposed, variance_g_proposed, p_x_proposed,
                                                       spike_allocation = inactive_spike_allocation[out_spike_idx])

    proposed_varianceGProbability <- .self$computeVarianceGPriorProbability(variance_g_proposed, Sg_proposed)

    Lp_pp <- proposed_XgLikelihood

    Lp_pp[inactive_outspike_idx] <- Lp_pp[inactive_outspike_idx] +
      .self$computeInactiveLikelihood(Yg_proposed[inactive_outspike_idx], spike_probability, inactive_means,
                                      inactive_variances, inactive_spike_allocation[out_spike_idx][inactive_outspike_idx])

    Lp_pp[active_outspike_idx] <- Lp_pp[active_outspike_idx] +
      .self$computeActiveLikelihood(Yg_proposed[active_outspike_idx], allocation_within_active[[1]][out_spike_idx][active_outspike_idx],
                                    active_means, active_variances)

    Lp_pp <- Lp_pp + proposed_varianceGProbability


    # accept or reject proposed move

    R = temperature * (Lp_pp - Lc_pc) + HR
    if(recover_x) recover()

    current_and_proposed_Yg <- cbind(Yg[out_spike_idx], Yg_proposed)
    current_and_proposed_variance_g <- cbind(variance_g[out_spike_idx], variance_g_proposed)

    current_and_proposed_likelihood <- cbind(XgLikelihood[out_spike_idx], proposed_XgLikelihood)


    choice_idx <- (log(runif(num_choice)) < R) + 1  # 1 means reject proposal, 2 means accept proposal.
    choice_matrix <- cbind(seq(num_choice), choice_idx) # row-column pairs for current_and_proposed_variance_g_probability with 1st column meaning reject proposal

    Yg[out_spike_idx] <<- current_and_proposed_Yg[cbind(seq(num_choice), choice_idx)]
    variance_g[out_spike_idx] <<- current_and_proposed_variance_g[choice_matrix]

    ## update XgLikelihood, variance_g, Sg and p_x

    XgLikelihood[out_spike_idx] <<- current_and_proposed_likelihood[choice_matrix]

    .self$set_sigmaX_pX()

    current_and_proposed_varianceGProbability <- cbind(variance_g_probability[out_spike_idx], proposed_varianceGProbability)

    variance_g_probability[out_spike_idx] <<- current_and_proposed_varianceGProbability[choice_matrix]

  }

)



