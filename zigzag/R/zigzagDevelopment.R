zigzag$methods(



  ######################################
  ####### MCMC under development #######
  ######################################

  mc3 = function(nchains = 4, burnin_ngen = 1000, mcmc_ngen = 10000, write_to_files = TRUE, mcmcprefix = "mc3_out", append = FALSE, recover_x = FALSE){

    if(append == FALSE){

      .self$mcmc(sample_frequency = 10, ngen = 1, write_to_files = TRUE, progress_plot = F, mcmcprefix = mcmcprefix, compute_probs = F, append = FALSE)

    }

    mcgen <- .self$gen

    proposal_probs[c(1,2,24)] <<- c(10, 0, 1)

    dT <- 0.025

    temps <- c()

    for(i in seq(nchains)) temps <- append(temps, 1/(1 + dT * (i - 1)))

    chain_list <- list()

    for(i in seq(nchains)){

      chain_list[[i]] <- .self$copy()

      chain_list[[i]]$temperature <- temps[i]

    }

    ## Burnin for all chains

    for(i in seq(nchains)){

      chain_list[[i]]$burnin(sample_frequency = 10, ngen = burnin_ngen, burnin_target_acceptance_rate=0.23, threads = 1, target_ESS = NULL,
                             progress_plot = FALSE, write_to_files = FALSE, burninprefix = paste0("burnin","_temp_",i))

      chain_list[[i]]$gen <- mcgen


    }



    ## run mcmc for x generations for each chain

    for(cycle in seq(mcmc_ngen)){

      for(i in seq(nchains)){

        if(chain_list[[i]]$temperature == 1){

          wtf <- TRUE

          apnd <- TRUE

        }else{

          wtf <- FALSE

          apnd <- FALSE

        }

        chain_list[[i]]$mcmc(sample_frequency = 20, ngen = 1, write_to_files = wtf, progress_plot = F, mcmcprefix = mcmcprefix, compute_probs = F, append = apnd)

        if(chain_list[[i]]$temperature == 1) chain_list[[i]]$plotDensities()

      }


      ## propose to swap each neighboring pair of chains

      # swapchains <- sample(seq(nchains - 1), size = 1, replace = FALSE)


      for(i in seq(nchains - 1)){

        #proposed chains to swap

        swapchains <- i

        posterior1_p <- chain_list[[swapchains]]$get_logPosterior(chain_list[[swapchains + 1]]$temperature)
        posterior2_p <- chain_list[[swapchains + 1]]$get_logPosterior(chain_list[[swapchains]]$temperature)

        posterior1_c <- chain_list[[swapchains]]$get_logPosterior()
        posterior2_c <- chain_list[[swapchains + 1]]$get_logPosterior()



        ## Make chain swap choice

        R = posterior1_p + posterior2_p - posterior1_c - posterior2_c


        if(log(runif(1)) < R){

          t1 <- chain_list[[swapchains]]$temperature
          t2 <- chain_list[[swapchains + 1]]$temperature

          ## swap temperatures, write_to_file boolean and allocation_trace
          chain_list[[swapchains]]$set_temperature(t2)
          chain_list[[swapchains + 1]]$set_temperature(t1)

          if(1 %in% c(t1, t2)){

            cold_idx <- c(swapchains, swapchains + 1)[which(c(t1, t2) == 1)]
            hot_idx <- c(swapchains, swapchains + 1)[which(c(t1, t2) != 1)]

            chain_list[[hot_idx]]$set_tuningParam_fields( chain_list[[cold_idx]]$get_tuningParam_fields() )

            chain_list[[hot_idx]]$allocation_trace <- chain_list[[cold_idx]]$allocation_trace


            # chain_list[[swapchains[hot]]]$wtf <<- chain_list[[swapchains[cold]]]$wtf # wtf is write to file boolean


          }

          chain_list[c(swapchains, swapchains+1)] <- chain_list[c(swapchains+1, swapchains)]

          cat(t1, "  ", t2, "\n")

          if(recover_x) recover()

        } #end proposal acceptence and parameter update

      } #end swap chain

    } #end mc3


  },


  estimateMarginalLikelihood = function(burnin_gen = 1000, mcmc_gen = 10000, ss_sample_frequency = 10, power_alpha = 0.2,
                                        stones = 10, cores = detectCores(), run_parallel=T, write_stones_to_files = F, ESS_target = 200){

    ssbeta <- ((1:stones)/stones)^(1/power_alpha)

    #### function to compute beta_k * lnl ####

    power_lnl <- function(k){

      dir.create(paste0(output_directory, "/", "marginalLikelihood"))

      require(MCMCpack)
      require(matrixStats)
      require(coda)

      beta <<- ssbeta[k]

      .self$burnin(sample_frequency = ss_sample_frequency, ngen = burnin_gen, write_to_files = F, target_ESS = floor(0.5*ESS_target))

      .self$mcmc(ngen = mcmc_gen, sample_frequency = ss_sample_frequency, write_to_files = write_stones_to_files,
                 mcmcprefix = paste0("marginalLikelihood/stone_", k), compute_probs = F, target_ESS = ESS_target)

      sslnl   = unlist(lnl_trace)/beta

      max_lnl = max(sslnl)

      ssloglike <- (ssbeta[k+1] - ssbeta[k]) * max_lnl + log(mean(exp((ssbeta[k+1] - ssbeta[k]) * (sslnl - max_lnl))))

      return(list(lnl_trace, beta, ssloglike))

    }

    ## Run stepping Stones.

    if(run_parallel){ # run parallel stepping stones

      cl <- makeCluster(cores)
      logmarginal <- clusterApply(cl, x=1:(length(ssbeta)-1), fun = power_lnl)

      stopCluster(cl)

    }else{ # run seriel stepping stones

      sslikelihoods <- list()

      for(k in 1:(length(ssbeta) - 1)){

        sslikelihoods[[length(sslikelihoods) + 1]] <- list(sslikelihoods, power_lnl(k)) #lnl_trace = beta * lnl

      }

      logmarginal <- sum(sslikelihoods[,2])

      beta <<- 1

    }

    sslnl_hat <- sapply(1:(stones - 1),function(s){return(logmarginal[[s]][[3]])})

    # blnl      <- as.matrix(sapply(1:(stones - 1),function(s){return(as.numeric(logmarginal[[s]][[1]]))}))
    # blnl      <- cbind(seq(nrow(blnl))-1,blnl)
    # write.table(round(blnl,6), file = paste0(output_directory, "/", "marginalLikelihood/beta_lnL.log"), quote = F, col.names = c("gen", paste0("stone",seq(length(ssbeta)-1),"_beta",round(ssbeta[-length(ssbeta)],5),"_lnL")), sep = "\t", row.names = F)

    return(sum(sslnl_hat))


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
    # p_x_proposed <- .self$xxget_px(yy = yg_proposed, gl = gene_lengths[zero_and_inactive_idx], xx = Xg[zero_and_inactive_idx,], sg = variance_g[zero_and_inactive_idx])


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
    # p_x_proposed <- .self$xxget_px(yy = Yg_proposed,  gl = gene_lengths[all_zero_idx], xx = Xg[all_zero_idx,], sg = variance_g[all_zero_idx])

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


  },

  ########################
  ## Plotting ############
  ########################

  plotPosteriorDensities = function(log_data_frame){

    # Lower level pdf file
    pdf("lower_level_posterior_densities.pdf")
    p_params = c("s0", "s1", "tau", "alpha_r")
    hyper_prior_list = list(c())
    layout = matrix(seq(4), ncol = 2)

    for(i in seq(4)){

      #plot posterior
      post = as.matrix(log_data_frame[, grepl(colnames(log_data_frame), pattern = p_params[i])])
      plot(density(post[,1]), lwd = 2, main = p_params[i])
      sapply(seq(ncol(post)), function(cc) lines(density(post[,cc])))

      #plot prior


    }
    dev.off()

    # Upper level pdf file
    pdf("Upper_level_posterior_densities.pdf")
    dev.off()

  },


  plotPosteriorTraces = function(log_data_frame){



  },



  ########################
  ## Depricated ##########
  ########################

  old_mhYg = function(recover_x = FALSE, tune = FALSE){

    Yg_proposed <<- Yg + rnorm(num_transcripts, 0, tuningParam_yg) * (1 - inactive_spike_allocation) #only change out of spike genes

    Sg_proposed <- .self$get_sigmax(yy = Yg_proposed)

    p_x_proposed <- .self$get_px(yy = Yg_proposed)

    ## proposed Yg likelihoods
    proposed_XgLikelihood = .self$computeXgLikelihood(Xg, Yg_proposed, variance_g, p_x_proposed)

    proposed_varianceGProbability = .self$computeVarianceGPriorProbability(variance_g, Sg_proposed)

    Lp_pp = proposed_XgLikelihood

    Lp_pp[inactive_idx] = Lp_pp[inactive_idx] +
      .self$computeInactiveLikelihood(Yg_proposed[inactive_idx], spike_probability, inactive_means,
                                      inactive_variances, inactive_spike_allocation[inactive_idx])

    Lp_pp[active_idx] = Lp_pp[active_idx] +
      .self$computeActiveLikelihood(Yg_proposed[active_idx],  allocation_within_active[[1]][active_idx],
                                    active_means, active_variances)

    Lp_pp = Lp_pp + proposed_varianceGProbability


    ## current Yg likelihoods
    Lc_pc = XgLikelihood

    Lc_pc[inactive_idx] = Lc_pc[inactive_idx] +
      .self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means,
                                      inactive_variances, inactive_spike_allocation[inactive_idx])

    Lc_pc[active_idx] = Lc_pc[active_idx] +
      .self$computeActiveLikelihood(Yg[active_idx],  allocation_within_active[[1]][active_idx],
                                    active_means, active_variances)

    Lc_pc = Lc_pc + variance_g_probability

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

    current_and_proposed_varianceGProbability <- cbind(variance_g_probability, proposed_varianceGProbability)

    variance_g_probability[out_spike_idx] <<- current_and_proposed_varianceGProbability[choice_matrix][out_spike_idx]

  }


)
