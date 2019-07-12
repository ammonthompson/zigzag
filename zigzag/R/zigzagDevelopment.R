zigzag$methods(



  ######################################3
  ####### under development #############
  #######################################

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



  }
)
