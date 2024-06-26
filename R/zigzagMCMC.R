zigzag$methods(

#' mcmc
#' @name mcmc
#' @description MCMC sampling from hieararchical bayesian mixture model
#' @usage my_zigzag_object$mcmc(sample_frequency = 100, ngen = 10000, target_ESS = NULL,
#' progress_plot = FALSE, write_to_files = TRUE, mcmcprefix = "out",
#' compute_probs = TRUE, run_posterior_predictive = FALSE, append = FALSE)
#' @param sample_frequency Number of generations between samples from the chain
#' @param ngen Number of generations to run the chain
#' @param target_ESS Run chain until the target effective sample size is achieved for the likelihood
#' @param write_to_files Write posterior samples to output files located in mcmcprefix_mcmc_ouput directory
#' @param mcmcprefix prefix for mcmc output directory as well as the output files within.
#' @param compute_probs compute probabilities of activity for each gene. Default is TRUE.
#' @param run_posterior_predictive run 500 posterior predictive simulation using parameters sampled from the chain.
#' @param run_posterior_predictive_and_plot run posterior predictive simulation and plot 50 simulated densities over observed densities.
#' @param append If posterior sample log files in mcmcprefix_mcmc_ouput already exists, then append samples to those files. Default is FALSE.
#'
#' @return Posterior log files including probability of active expression and MCMC diagnostics in mcmc_report directory.
#' @export
#'
#' @examples{
#' \dontrun{mcmc(sample_frequency = 100, ngen = 50000, mcmcprefix = "organ_x_posterior_sample")}
#' }
#'
  mcmc = function(sample_frequency                  = 100,
                  ngen                              = 10000,
                  target_ESS                        = NULL,
                  write_to_files                    = TRUE,
                  mcmcprefix                        = "out",
                  compute_probs                     = TRUE,
                  run_posterior_predictive          = FALSE,
                  run_posterior_predictive_and_plot = run_posterior_predictive,
                  append                            = FALSE){


    ###########################
    # initialize output files #
    ###########################
    mcmc_prefixdir <- paste0(mcmcprefix, "_mcmc_output")
    vprefix <- unlist(strsplit(mcmcprefix, split = "/"))
    prefix <- vprefix[length(vprefix)]


    if(write_to_files & append == FALSE){

      dir.create(paste0(output_directory, "/", mcmc_prefixdir))
      .self$initializeOutputFiles(paste0(mcmc_prefixdir, "/", prefix),
                                  (run_posterior_predictive | run_posterior_predictive_and_plot))


    }else if(append & !file.exists(paste0(output_directory, "/",
                                          mcmc_prefixdir, "/",
                                          prefix, "_model_parameters.log"))){

      print("No file exists to append. Creating new files")
      print(paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, "_model_parameters.log"))
      dir.create(paste0(output_directory, "/", mcmc_prefixdir))
      .self$initializeOutputFiles(paste0(mcmc_prefixdir, "/", prefix),
                                  (run_posterior_predictive | run_posterior_predictive_and_plot))

    }

    # if(length(dev.list()) == 0 & run_posterior_predictive_and_plot) dev.new()


    #########################################################
    # Initialize posterior predictive density plot of       #
    # Xg level 1 (lower level) and Yg level 2 (upper level) #
    #########################################################
    posXg_rowMeans <- apply(Xg, 1, function(g) if(max(g) > -Inf) mean(g[g > -Inf]) else -Inf)
    xg_colmeans <- apply(Xg, 2, function(x) mean(x[x > -Inf]))
    grandmean_xg_colmeans <- mean(xg_colmeans)
    scaled_xg <- sapply(seq(num_libraries), function(lib) Xg[,lib] - xg_colmeans[lib] + grandmean_xg_colmeans)
    num_libs_plot_postpred <- num_libraries
    if(num_libs_plot_postpred > 25) num_libs_plot_postpred = 25

    plot_rows <- sqrt(num_libs_plot_postpred) - sqrt(num_libs_plot_postpred)%%1 + (sqrt(num_libs_plot_postpred)%%1 > 0)
    plot_cols <- num_libs_plot_postpred/plot_rows - (num_libs_plot_postpred/plot_rows)%%1 + ((num_libs_plot_postpred/plot_rows)%%1 > 0)
    multi_plot_pars <- list()

    # Set a maximum of 500 posterior predictive simulations, and 50 post. pred. density plots
    postpred_freq <- round(ngen/(500 * sample_frequency))
    if(postpred_freq == 0) postpred_freq <- 1
    postpred_plotfeq <- round(ngen/(50 * postpred_freq * sample_frequency))
    if(postpred_plotfeq == 0) postpred_plotfeq <- 1

    if(run_posterior_predictive_and_plot){

      ## all libs up to 25 for lower level
      pdf(file = paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, ".scaled_L1_post_pred.pdf"))
      par("mfrow" = c(plot_rows, plot_cols))
      par("mai" = c(0.4, 0.2, 0.2, 0.2))

      for(plib in seq(num_libs_plot_postpred)){

        libdensity = density(scaled_xg[, plib])
        plot(libdensity$x, libdensity$y, type = "l", ylim = c(0, max(libdensity$y)*1.1),
             lwd = 2, ylab = "", xlab = "", main = paste0(mcmcprefix, " Lib. ", plib), axes = F)
        axis(1)
        multi_plot_pars[[plib]] = par("usr","plt","mfg")

      }

      post_pred_multi_L1_plot_device = dev.cur()


      ## Upper level plot
      d_posXg_rowMeans = density(posXg_rowMeans)
      d_posYg_rowMeans = density(Yg[out_spike_idx])

      post_pred_plot_xlims = c(min(Yg[out_spike_idx]) - 1, max(Yg[out_spike_idx]) + 1)

      pdf(file = paste0(output_directory, "/", mcmc_prefixdir, "/", prefix, ".L2_post_pred.pdf"))

      y_max =  max(d_posYg_rowMeans$y)

      layout(matrix(seq(2),nrow=2))
      plot(NULL, xlim = post_pred_plot_xlims, main = paste0(mcmcprefix, ", true expression (Y) distribution"),
           ylim = c(0, 1.2 * y_max),
           xlab = "log Expression", ylab = "density")

      legend(post_pred_plot_xlims[1], 1.2 * y_max, legend = c("Post. Y", "Sim. Y"),
             col = c("orange", "green"), lty = 1, lwd = 3)
      post_pred_L2_plot_device = dev.cur()

      plot(NULL, xlim = post_pred_plot_xlims, main = "",
           ylim = c(0, 1.2 * y_max),
           xlab = "log Expression", ylab = "density")
      legend(post_pred_plot_xlims[1], 1.2 * y_max, legend = c("Post. inactive Y", "Post. active Y"),
             col = c("dodgerblue", "tomato"), lty = 1, lwd = 3)
      post_L2_plot_device = dev.cur()



    } # end post predictive plotting setup


    #########################
    ### Set up MCMC params###
    #########################
    all_allocation_mcmc <- allocation_active_inactive * allocation_within_active[[1]]

    .self$allocation_trace <- matrix(apply(component_matrix, 2, function(comp_matrix_col) 1 *
                        (comp_matrix_col == all_allocation_mcmc)), nrow = num_transcripts)

    .self$lnl_trace[[length(lnl_trace) + 1]] <- .self$calculate_lnl(num_libraries)

    # j increments with gen when target_ESS isn't used
    # i <- gen
    # j <- 0
    j <- gen

    plist_length <- length(proposal_list)

    yg_varg_subsample_frequency <- floor(ngen/(sample_frequency * 500))
    if(yg_varg_subsample_frequency == 0) yg_varg_subsample_frequency <- 1

    if(write_to_files) write_yg_varg = TRUE

    starting_gen <- gen
    ngen <- ngen + gen


    ##############
    # begin mcmc #
    ##############

    while(j <= ngen ){

      ### run n = plist_length proposals ###
      for(p in sample(seq(plist_length), plist_length, replace = TRUE, prob = proposal_probs))
        proposal_list[[p]]()

      #############################
      # sample from markov chain  #
      #############################
      # if(i %% sample_frequency == 0){
      #
      #   if(is.null(target_ESS))
      #     j <- i

        if(is.null(target_ESS))
          j <- gen

        gen_write <- FALSE

        if(gen %% sample_frequency == 0){

          # cat("#### ",i," ####", "\n")
          cat("#### ",gen," ####", "\n")

        ### For computing Pr(zag = 1) ###
        if(temperature == 1){

          all_allocation_mcmc <- allocation_active_inactive * allocation_within_active[[1]]
          .self$allocation_trace <- allocation_trace + apply(component_matrix, 2, function(comp_matrix_col) 1 *
                                                          (comp_matrix_col == all_allocation_mcmc))

        }

        .self$lnl_trace[[length(lnl_trace) + 1]] <- .self$calculate_lnl(num_libraries)


        if(write_to_files & write_yg_varg * temperature == 1){

          # .self$writeToOutputFiles(paste0(mcmc_prefixdir,"/", mcmcprefix), gen = i)
          .self$writeToOutputFiles(paste0(mcmc_prefixdir,"/", mcmcprefix), gen = gen)

          # if((i / sample_frequency) %% yg_varg_subsample_frequency == 0)
          #     .self$writeToYgVariancegOutputFiles(paste0(mcmc_prefixdir,"/", mcmcprefix), gen = i)
          # if((i - starting_gen) / (yg_varg_subsample_frequency * sample_frequency) > 500)
          #   write_yg_varg <- FALSE
          if((gen / sample_frequency) %% yg_varg_subsample_frequency == 0)
            .self$writeToYgVariancegOutputFiles(paste0(mcmc_prefixdir,"/", mcmcprefix), gen = gen)
          if((gen - starting_gen) / (yg_varg_subsample_frequency * sample_frequency) > 500)
            write_yg_varg <- FALSE

          gen_write <- TRUE
        }


        if(!is.null(target_ESS) & length(lnl_trace) > 100){

          if(.self$calculate_lnl_ESS() > target_ESS) break

        }

      } #end sample from chain



      #######################################
      # Run posterior predictive simulation #
      #######################################

    # if((run_posterior_predictive | run_posterior_predictive_and_plot) & (i / (postpred_freq * sample_frequency)) %% 1 == 0){
    #
    #   post_pred_instance <- .self$posteriorPredictiveSimulation(prefix = mcmcprefix)
    #
    #   if( run_posterior_predictive_and_plot & (i / (postpred_freq * postpred_plotfeq * sample_frequency)) %% 1 == 0){
    #
    #     .self$update_postPredPlots(post_pred_instance, multi_plot_pars,
    #                                post_pred_multi_L1_plot_device,
    #                                post_pred_L2_plot_device, d_posXg_rowMeans)
    #
    #   } #end posterior predicitve plotting
    #
    # } #end posterior predictive simulation
    if((run_posterior_predictive | run_posterior_predictive_and_plot) & (gen / (postpred_freq * sample_frequency)) %% 1 == 0){

      post_pred_instance <- .self$posteriorPredictiveSimulation(prefix = mcmcprefix)

      if( run_posterior_predictive_and_plot & (gen / (postpred_freq * postpred_plotfeq * sample_frequency)) %% 1 == 0){

        .self$update_postPredPlots(post_pred_instance, multi_plot_pars,
                                   post_pred_multi_L1_plot_device,
                                   post_pred_L2_plot_device, d_posXg_rowMeans)

      } #end posterior predicitve plotting

    } #end posterior predictive simulation



    # i <- i + 1
    # .self$gen <- i
    .self$gen <- gen + 1

      #if postpred interrupted, close out pdf devices
      on.exit(if(length(dev.list()) > 2) dev.off(post_pred_multi_L1_plot_device))
      on.exit(if(length(dev.list()) > 2) dev.off(post_pred_L2_plot_device), add = TRUE)
      on.exit(if(gen_write) .self$gen <- gen + 1, add = TRUE, after = TRUE)

    } #end mcmc


    ############################
    # finalize post pred plots #
    ############################
    if(run_posterior_predictive_and_plot){

      .self$finalize_postPredPlots(multi_plot_pars,
                                   post_pred_multi_L1_plot_device,
                                   post_pred_L2_plot_device,
                                   d_posXg_rowMeans, scaled_xg)
    }

    #####################################
    # compute gene activity state      ##
    # probabilities and write to files ##
    #####################################
    if(compute_probs && write_to_files)
      .self$computeGeneExpressionProbs_writeToFile(paste0(mcmc_prefixdir,"/", prefix))


    #####################################
    # create MCMC diagnostic reports   ##
    #####################################
    .self$generate_MCMC_reports(paste0(output_directory, "/", mcmc_prefixdir), prefix)


    on.exit(if(length(dev.list()) > 2) dev.off(post_pred_multi_L1_plot_device))
    on.exit(if(length(dev.list()) > 2) dev.off(post_pred_L2_plot_device), add = TRUE)

    # get the rest of em!
    while(dev.cur() > 1) {
      dev.off()
    }

  }

)
