zigzag$methods(

#' mcmc
#' @name mcmc
#' @description MCMC sampling from hieararchical bayesian mixture model
#' @usage my_zigzag_object$mcmc(sample_frequency = 10, ngen = 10000, target_ESS = NULL,
#' progress_plot = FALSE, write_to_files = TRUE, mcmcprefix = "out",
#' compute_probs = TRUE, run_posterior_predictive = TRUE, append = FALSE)
#' @param sample_frequency Number of generations between samples from the chain
#' @param ngen Number of generations to run the chain
#' @param target_ESS Run chain until the target effective sample size is achieved for the likelihood
#' @param progress_plot Display a plot showing summary statistics and graphs. (in Rstudio)
#' @param write_to_files Write posterior samples to output files located in mcmcprefix_mcmc_ouput directory
#' @param mcmcprefix prefix for mcmc output directory as well as the output files within.
#' @param compute_probs compute probabilities of activity for each gene. Default is TRUE.
#' @param run_posterior_predictive run posterior predictive simulation for every other sample from the chain.
#' @param append If posterior sample log files in mcmcprefix_mcmc_ouput already exists, then append samples to those files. Default is FALSE.
#'
#' @return blorb
#' @export
#'
#' @examples{
#' \dontrun{mcmc(sample_frequency = 100, ngen = 50000, mcmcprefix = "organ_x_posterior_sample")}
#' }
#'
  mcmc = function(sample_frequency = 10, ngen = 10000, target_ESS = NULL, progress_plot = FALSE,
                  write_to_files = TRUE, mcmcprefix = "out", compute_probs = TRUE,
                  run_posterior_predictive = TRUE, append = FALSE){

    ###########################
    # initialize output files #
    ###########################
    mcmc_prefixdir = paste0(mcmcprefix, "_mcmc_output")
    vprefix = unlist(strsplit(mcmcprefix, split = "/"))
    prefix = vprefix[length(vprefix)]


    if(write_to_files & append == FALSE){

      dir.create(paste0(output_directory, "/", mcmc_prefixdir))
      .self$initializeOutputFiles(paste0(mcmc_prefixdir, "/", prefix), run_posterior_predictive)


    }else if(append & !file.exists(paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, "_model_parameters.log"))){

      print("No file exists to append. Creating new files")
      print(paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, "_model_parameters.log"))
      dir.create(paste0(output_directory, "/", mcmc_prefixdir))
      .self$initializeOutputFiles(paste0(mcmc_prefixdir, "/", prefix), run_posterior_predictive)

    }


    #############################################
    # Set up progress plotting on active device #
    #############################################
    if(progress_plot){

      if(length(dev.list()) == 0) dev.new()
      progress_plot_device = dev.list()[1]

      histogram_list = list()
      histogram_list[[1]] <- hist(Yg, breaks = 100, plot = FALSE)

      max_heights = max(histogram_list[[1]]$density[-1])
      total_height = sum(max_heights) - 0.025
      max_height = max(max_heights)

      baselines = c(0)

      plot(NULL, xlim=c(-10,10), ylim = c(0, 1.01 * total_height))
      .self$drawLibraryHistograms(histogram_list[[1]]$mids, histogram_list[[1]]$density, bl = baselines[1])

    }

    #########################################################
    # Begin posterior predictive density plot of            #
    # Xg level 1 (lower level) and Yg level 2 (upper level) #
    #########################################################
    posXg_rowMeans = apply(Xg, 1, function(g) if(max(g) > -Inf) mean(g[g > -Inf]) else -Inf)
    xg_colmeans = apply(Xg, 2, function(x) mean(x[x > -Inf]))
    grandmean_xg_colmeans = mean(xg_colmeans)
    scaled_xg = sapply(seq(num_libraries), function(lib) Xg[,lib] - xg_colmeans[lib] + grandmean_xg_colmeans)

    num_libs_plot_postpred = num_libraries
    if(num_libs_plot_postpred > 25) num_libs_plot_postpred = 25

    plot_rows = sqrt(num_libs_plot_postpred) - sqrt(num_libs_plot_postpred)%%1 + (sqrt(num_libs_plot_postpred)%%1 > 0)
    plot_cols = num_libs_plot_postpred/plot_rows - (num_libs_plot_postpred/plot_rows)%%1 + ((num_libs_plot_postpred/plot_rows)%%1 > 0)

    multi_plot_pars = list()

    if(run_posterior_predictive){

      ## all libs up to 25 for lower level
      pdf(file = paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, ".scaled_L1_post_pred.pdf"))

      par("mfrow" = c(plot_rows, plot_cols))

      par("mai" = c(0.4, 0.2, 0.2, 0.2))

      for(plib in seq(num_libs_plot_postpred)){
        libdensity = density(scaled_xg[, plib])
        plot(libdensity$x, libdensity$y, type = "l", ylim = c(0, max(libdensity$y)*1.1),
             lwd = 2, ylab = "", xlab = "", main = paste0("Lib. ", plib), axes = F)
        axis(1)
        multi_plot_pars[[plib]] = par("usr","plt","mfg")
      }

      post_pred_multi_L1_plot_device = dev.cur()


      ## Upper level plot
      d_posXg_rowMeans = density(posXg_rowMeans)
      d_posYg_rowMeans = density(Yg[out_spike_idx])

      post_pred_plot_xlims = c(min(posXg_rowMeans[posXg_rowMeans > -Inf]) - 1, max(posXg_rowMeans) + 1)

      pdf(file = paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, ".L2_post_pred.pdf"))

      y_max =  max(c(d_posXg_rowMeans$y, d_posYg_rowMeans$y))

      plot(NULL, xlim = post_pred_plot_xlims, main = "Upper level",
           ylim = c(-1.2 * y_max, 1.2 * y_max),
           xlab = "log Expression", ylab = "density")
      abline(h=0, lwd = 2)

      legend(post_pred_plot_xlims[1], 1.2 * y_max, legend = c("Post. Y", "Sim. Y"),
             col = c("orange", "green"), lty = 1, lwd = 3)
      legend(post_pred_plot_xlims[1],  -y_max, legend = c("Post.Y - Sim.Y"),
             col = c("red"), lty = 1, lwd = 3)

      post_pred_L2_plot_device = dev.cur()

     }

    if(progress_plot) dev.set(progress_plot_device)


    ### Set up MCMC ###
    all_allocation_mcmc = allocation_active_inactive * allocation_within_active[[1]]

    allocation_trace <<- matrix(apply(component_matrix, 2, function(comp_matrix_col) 1 *
                        (comp_matrix_col == all_allocation_mcmc)), nrow = num_transcripts)

    lnl_trace[[length(lnl_trace) + 1]] <<- .self$calculate_lnl(num_libraries)

    i <- gen
    j <- 0

    proposal_list <- list(.self$gibbsMixtureWeights,
                          .self$gibbsAllocationActiveInactive,
                          .self$gibbsAllocationWithinActive,
                          .self$mhInactiveMeans,
                          .self$mhInactiveVariances,
                          .self$mhActiveMeansDif,
                          .self$mhActiveVariances,
                          .self$gibbsSpikeProb,
                          .self$mhSpikeAllocation,
                          .self$mhYg,
                          .self$mhSigma_g,
                          .self$mhTau,
                          .self$mhSg,
                          .self$mhS0Tau,
                          .self$mhP_x )


    plist_length <- length(proposal_list)

    meanhistory = c()

    ##############
    # begin mcmc #
    ##############

    while(j <= ngen ){

      ### run n = plist_length proposals ###
      for(p in sample(seq(plist_length), plist_length, replace = TRUE, prob = proposal_probs)) proposal_list[[p]]()

      #############################
      # sample from markov chain  #
      #############################
      if(i %% sample_frequency == 0){

        if(is.null(target_ESS)) j = i

        cat("#### ",i," ####", "\n")

        if(progress_plot) .self$drawLibraryDensities(histogram_list[[1]]$mids, 1, baselines[1], alpha = 0.05)

        ### For computing Pr(zag = 1) ###
        if(temperature == 1){

          all_allocation_mcmc = allocation_active_inactive * allocation_within_active[[1]]
          allocation_trace <<- allocation_trace + apply(component_matrix, 2, function(comp_matrix_col) 1 * (comp_matrix_col == all_allocation_mcmc))

        }

        lnl_trace[[length(lnl_trace) + 1]] <<- .self$calculate_lnl(num_libraries)


        if(write_to_files & temperature == 1) .self$writeToOutputFiles(paste0(mcmc_prefixdir,"/", mcmcprefix), gen = i)

        if(!is.null(target_ESS) & length(lnl_trace) > 100){

          if(.self$calculate_lnl_ESS() > target_ESS) break

        }

        #######################################
        # Run posterior predictive simulation #
        #######################################
        #add a red density line to the posterior predictive plot for simulated mean data|model params
        if(run_posterior_predictive & (i / (2 * sample_frequency)) %% 1 == 0){

          post_pred_instance = .self$posteriorPredictiveSimulation(prefix = mcmcprefix)

          if((i / (2 * 2 * sample_frequency)) %% 1 == 0){
            cat("plotting post pred...\n")

            ### Multiple libs plots ###
            dev.set(post_pred_multi_L1_plot_device)

            lib_means_sim_xg = apply(post_pred_instance[[1]], 2, function(sx) mean(sx[sx > -Inf]))
            grand_mean_sim_xg = mean(lib_means_sim_xg)
            scaled_sim_xg = sapply(seq(num_libraries), function(sx) post_pred_instance[[1]][,sx] - lib_means_sim_xg[sx] + grand_mean_sim_xg)

            for(lib_plot in seq(num_libs_plot_postpred)){
              par("usr" = multi_plot_pars[[lib_plot]]$usr)
              par("plt" = multi_plot_pars[[lib_plot]]$plt)
              par("mfg" = multi_plot_pars[[lib_plot]]$mfg)
              lines(density(scaled_sim_xg[,lib_plot]), col = rgb(1,0,0,0.1))
            }

            ### Upper level plot ###
            dev.set(post_pred_L2_plot_device)

            ybins = seq(min(c(Yg[out_spike_idx] , post_pred_instance[[2]][out_spike_idx])) -1 ,
                        max(c(Yg[out_spike_idx] , post_pred_instance[[2]][out_spike_idx]))+1, by = 0.2)

            dYg = hist(Yg[out_spike_idx], plot = F, breaks = ybins)
            dsimYg = hist(post_pred_instance[[2]][out_spike_idx], plot = F, breaks = ybins)

            lines(lowess(dYg$mids, dYg$density - dsimYg$density - 0.6 * max(d_posXg_rowMeans$y), f=0.1) , col=rgb(1,0,0,0.05))
            lines(density(Yg[out_spike_idx]), col = rgb(1,0.5,0,0.05))
            lines(density(post_pred_instance[[2]][out_spike_idx]), col = rgb(0,0.85,0,0.05))

            if(progress_plot) dev.set(progress_plot_device)

          } #end posterior predicitve simulation

        } #end posterior predictive processing

      } #end sample from chain

      i = i + 1

      gen <<- i

      on.exit(if(length(dev.list()) > 2) dev.off(post_pred_multi_L1_plot_device))
      on.exit(if(length(dev.list()) > 2) dev.off(post_pred_L2_plot_device), add = TRUE)

    } #end mcmc


    ############################
    # finalize post pred plots #
    ############################
    if(run_posterior_predictive){
      #Level 1
      dev.set(post_pred_multi_L1_plot_device)
      for(i in seq(num_libs_plot_postpred)){
        par("usr" = multi_plot_pars[[i]]$usr)
        par("plt" = multi_plot_pars[[i]]$plt)
        par("mfg" = multi_plot_pars[[i]]$mfg)
        lines(density(scaled_xg[,i]), lwd = 0.75)
      }

      #Upper level (2)
      dev.set(post_pred_L2_plot_device)
      abline(h= -0.6 * max(d_posXg_rowMeans$y), lwd = 1.5)


      dev.off(post_pred_multi_L1_plot_device)
      dev.off(post_pred_L2_plot_device)
    }

    #####################################
    # compute gene activity state      ##
    # probabilities and write to files ##
    #####################################
    if(compute_probs && write_to_files) .self$computeGeneExpressionProbs_writeToFile(paste0(mcmc_prefixdir,"/", prefix))

    on.exit(if(length(dev.list()) > 2) dev.off(post_pred_multi_L1_plot_device))
    on.exit(if(length(dev.list()) > 2) dev.off(post_pred_L2_plot_device), add = TRUE)

  },



  #######################################################
  # Posterior predictive simulation to be called during #
  # the mcmc when run_posterior_predictive == TRUE      #
  #######################################################

  posteriorPredictiveSimulation = function(num_sims = 1000, plotlib = 1, prefix = "out", recover_x = FALSE){

    pp_prefix = paste0(output_directory, "/", prefix, "_mcmc_output/", prefix)

    ##################
    ### FUNCTIONS ####
    ##################

    #function pars
    binwidth = 0.1
    model_range = c(qnorm(0.000001, inactive_means, sqrt(inactive_variances)),
                    qnorm(0.999999, active_means[num_active_components], sqrt(active_variances[num_active_components])))
    Yg_range = c(min(Yg[out_spike_idx]), max(Yg[out_spike_idx]))
    xg_greater_minusInf_idx = which(Xg > -Inf, arr.ind = T)
    Xg_range = c(min(Xg[xg_greater_minusInf_idx]), max(Xg[xg_greater_minusInf_idx]))
    all_ranges = c(model_range, Yg_range, Xg_range)

    bins = seq(min(all_ranges) - 2, max(all_ranges) + 2, by = binwidth)
    yg_lessThan_bin = sapply(bins, function(y) (Yg < y) * 1 )

    get_cumulative = function(xdat){ #xdat is a vector (one library)

      zeroIdx = which(xdat == -Inf)
      num_zeros = length(zeroIdx)
      p_zeros = num_zeros/num_transcripts

      xdat = xdat[-zeroIdx]

      F = c()
      for( x in bins ) F = c( F, count( xdat <= x )/(num_transcripts - num_zeros))

      return( c(p_zeros, F) )

    }

    get_rumsfeld = function(xdat, rec_x = FALSE){

      model_Br = sapply(seq(num_libraries), function(r) (1 - spike_probability * (1 - weight_active)) * mean(1 - p_x[out_spike_idx, r]) +
                          spike_probability * (1 - weight_active))

      xdat_Br = sapply(seq(num_libraries), function(r) sum(xdat[,r] == -Inf)/num_transcripts)

      if(rec_x == T) recover()

      return(mean(abs(xdat_Br - model_Br)))

    }

    get_level1_wasserstein = function(xdat, rec_x = FALSE){

      model_cum = sapply(seq(num_libraries), function(r){
        colMeans(p_x[out_spike_idx, r] * yg_lessThan_bin[out_spike_idx, ])
      })


      norm_factor = colMaxs(model_cum)

      model_cum = t(t(model_cum)/norm_factor)

      xdat_cum = sapply(seq(num_libraries), function(r){
        sapply(bins, function(y) mean(xdat[xdat[,r] > -Inf, r] < y))
      })

      #average number of genes detected < x averaged over all libraries and the wasserstein distance.
      mean_xdat_cum = rowMeans(xdat_cum)
      mean_model_cum = rowMeans(model_cum)
      Wass_means = trapz_integration(bins, abs(mean_xdat_cum - mean_model_cum))

      if(rec_x == T) recover()

      return(c(Wass_means))

    }

    get_level2_wasserstein = function(ydat, rec_x = FALSE){

      model_cum = get_model_cumulative(bins)
      ydat_cum = sapply(bins[-1], function(x) mean(ydat[out_spike_idx] < x))

      Wass_means = trapz_integration(bins[-1], abs(ydat_cum - model_cum))

      if(rec_x == T) recover()

      return(Wass_means)

    }


    #####################
    ### END FUNCTIONS ###
    #####################


    #### Simulate data from posterior and compute realized discrepancy statstics

    sim_xg = sapply(1:num_libraries, function(x){return(rnorm(length(Yg), Yg, sqrt(sigma_g)))})

    # p detect
    sim_xg = sapply(1:num_libraries, function(lib){

      temp_sim_xg = sim_xg[,lib]
      not_detected = (runif(length(Yg)) < (1 - p_x[,lib]))
      temp_sim_xg[not_detected] <- -Inf

      return(temp_sim_xg)

    })

    ### simulate spike p( spike | inactive )
    sim_xg[in_spike_idx,] = rep(-Inf, num_libraries)

    sim_xg = as.data.frame(sim_xg)

    # level 1 discrepencies, zeros and non-zero Xg distros
    W_L1 = get_level1_wasserstein(Xg) - get_level1_wasserstein(sim_xg)
    B_rumsfeld = get_rumsfeld(Xg) - get_rumsfeld(sim_xg)


    ### level 2 discrepency simulating Yg
    #simulate Yg
    sim_yg = NULL
    sim_yg[inactive_idx] = rnorm(length(inactive_idx), inactive_means, sqrt(inactive_variances))
    sim_yg[active_idx] = rnorm(length(active_idx), active_means[allocation_within_active[[1]][active_idx]],
                               sqrt(active_variances[allocation_within_active[[1]][active_idx]]))


    W_L2 = get_level2_wasserstein(Yg) - get_level2_wasserstein(sim_yg)

    cat("W_L1: ",W_L1, "  W_L2: ", W_L2, '\n')

    if(recover_x) recover()

    write.table(matrix(c(gen, W_L1, W_L2, B_rumsfeld), nrow = 1),
                file = paste0(pp_prefix, ".post_predictive_output.log"),
                append=T, sep="\t", row.names=F, col.names=F)

    return(list(sim_xg, sim_yg))

  }


)
