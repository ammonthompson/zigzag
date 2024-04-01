zigzag$methods(

#' burnin
#' @name burnin
#' @description Run burnin and tune proposal size parameters for hieararchical bayesian mixture model
#' @usage burnin(sample_frequency = 100, ngen = 1000, burnin_target_acceptance_rate=0.44,
#' progress_plot = FALSE, write_to_files = TRUE, burninprefix = "burnin", append = FALSE)
#' @param sample_frequency Number of generations between samples from the chain
#' @param ngen Number of generations to run the chain
#' @param burnin_target_acceptance_rate proportion of proposals that are accepted. For proposal size tuning.
#' @param progress_plot Show plots of model statistics as MCMC progresses
#' @param write_to_files Write burnin samples to output files in burninprefix_burnin_output directory
#' @param burninprefix The prefix for the burnin output directory as well as the burnin output files within
#' @param append If files already exist from previous interrupted burnin, append to those files if TRUE.
#'
# @return
# @export
#'
# @examples
#'
  burnin = function(sample_frequency = 100,
                    ngen = 10000,
                    burnin_target_acceptance_rate=0.44,
                    progress_plot = FALSE,
                    write_to_files = TRUE,
                    burninprefix = "output",
                    append = FALSE){

  timestart = as.numeric(Sys.time())

  ###########################
  # initialize output files #
  ###########################
  mcmc_prefixdir = paste0(burninprefix, "_burnin")
  vprefix = unlist(strsplit(burninprefix, split = "/"))
  prefix = vprefix[length(vprefix)]

  if(write_to_files & append == FALSE){

   dir.create(paste0(output_directory, "/", mcmc_prefixdir))
   .self$initializeOutputFiles(paste0(mcmc_prefixdir, "/", prefix))

  }else if(append == TRUE & !file.exists(paste0(output_directory, "/",
                                                mcmc_prefixdir, "/",
                                                burninprefix,
                                                "_burnin_model_parameters.log"))){

   print("File does not exist to append. Creating new file")
   dir.create(paste0(output_directory, "/", mcmc_prefixdir))
   .self$initializeOutputFiles(paste0(mcmc_prefixdir, "/", prefix))

  }


  i = gen
  j = 0

  plist_length <- length(proposal_list)

  ################
  # run the mcmc #
  ################

  while(j <= ngen ){
   end = FALSE
   for(p in sample(seq(plist_length), plist_length, replace = TRUE, prob = proposal_probs)){
     proposal_list[[p]](tune = TRUE);

     if(is.na(sum(Yg)) || is.na(sum(variance_g)) ||is.na(sum(tuningParam_yg)) || is.na(sum(inactive_spike_allocation)) ||
        is.nan(sum(inactive_spike_allocation)) || is.na(sum(allocation_active_inactive))||
        is.nan(sum(allocation_active_inactive)) || is.na(sum(variance_g_trace)) || is.na(sum(Yg_trace))){

       print("an NA or NaN somewhere!")
       print(p)
       end = TRUE
       break

     }
   }
   if(end) break

   j <- i

   #####################
   # sample from chain #
   #####################
   if(i %% sample_frequency == 0 & i > 2 * sample_frequency){

     xx_timestart = as.numeric(Sys.time())


     cat("#### ",i," ####  ", lnl_trace[[length(lnl_trace)]], "\n")

     if(progress_plot) .self$burninProgressPlot()


     #####################################################################
     ## set trace vars (this doesn't store the values, just updates it. ##
     ## Important in case user interrupts the burnin)                   ##
     #####################################################################

     .self$s0_trace[[2]] <- s0
     .self$s1_trace[[2]] <- s1
     .self$alpha_r_trace[[2]] <- alpha_r

     .self$inactive_means_trace[[2]] <- inactive_means
     .self$inactive_variances_trace[[2]] <- inactive_variances
     .self$spike_probability_trace[[2]] <- spike_probability
     .self$active_means_trace[[2]] <- active_means
     .self$active_variances_trace[[2]] <- active_variances
     .self$weight_active_trace[[2]] <- weight_active
     .self$weight_within_active_trace[[2]] <- weight_within_active

     if(temperature == 1){

       all_allocation_burnin = allocation_active_inactive * allocation_within_active[[1]]

       .self$allocation_trace <- matrix(apply(component_matrix, 2, function(comp_matrix_col) 1 *
                                           (comp_matrix_col == all_allocation_burnin)), nrow = num_transcripts)

     }

     .self$lnl_trace[[length(lnl_trace) + 1]] <- .self$calculate_lnl(num_libraries)

     if(write_to_files & temperature == 1){

       .self$writeToOutputFiles(paste0(mcmc_prefixdir,"/", prefix), gen = i)

       if((i / sample_frequency) %% 4 == 0)
         .self$writeToYgVariancegOutputFiles(paste0(mcmc_prefixdir,"/", prefix), gen = i)

     }

     cat("time: ", as.numeric(Sys.time()) - timestart, "\n")
     timestart = as.numeric(Sys.time())

   }

   i = i + 1

   .self$gen <- i

   ################################################
   ### tune the tuning parameters (tuningParam) ###
   ################################################
   if( gen %% 400 == 0 ) .self$tune_all(burnin_target_acceptance_rate)

  } #end mcmc

  .self$lnl_trace <- list(lnl_trace[[length(lnl_trace)]])

  }


)
