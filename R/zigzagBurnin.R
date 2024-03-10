zigzag$methods(

#' burnin
#' @name burnin
#' @description Run burnin and tune proposal size parameters for hieararchical bayesian mixture model
#' @usage burnin(sample_frequency = 100, ngen = 1000, burnin_target_acceptance_rate=0.44,
#' target_ESS = NULL, progress_plot = FALSE, write_to_files = TRUE,
#' burninprefix = "burnin", append = FALSE)
#' @param sample_frequency Number of generations between samples from the chain
#' @param ngen Number of generations to run the chain
#' @param burnin_target_acceptance_rate proportion of proposals that are accepted. For proposal size tuning.
#' @param threads Depricated
#' @param target_ESS Depricated. Probably
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
  burnin = function(sample_frequency = 100, ngen = 10000, burnin_target_acceptance_rate=0.44,
                 threads = 1, target_ESS = NULL, progress_plot = FALSE, write_to_files = TRUE,
                 burninprefix = "output", append = FALSE){

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

  }else if(append == TRUE & !file.exists(paste0(output_directory, "/", mcmc_prefixdir, "/", burninprefix, "_burnin_model_parameters.log"))){

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

     if(min(XgLikelihood) == -Inf){print(p); recover()}

     if(is.na(sum(Yg)) || is.na(sum(variance_g)) ||is.na(sum(tuningParam_yg)) || is.na(sum(inactive_spike_allocation)) ||
        is.nan(sum(inactive_spike_allocation)) || is.na(sum(allocation_active_inactive))||
        is.nan(sum(allocation_active_inactive)) || is.na(sum(variance_g_trace[[1]])) || is.na(sum(Yg_trace[[1]]))){

       print("an NA or NaN somewhere!")
       print(p)
       end = TRUE
       break

     }
   }
   if(end) break


   #####################
   # sample from chain #
   #####################
   if(i %% sample_frequency == 0 & i > 2 * sample_frequency){

     xx_timestart = as.numeric(Sys.time())

     if(is.null(target_ESS)) j = i

     cat("#### ",i," ####  ", lnl_trace[[length(lnl_trace)]], "\n")

     if(progress_plot) .self$burninProgressPlot()


     #####################################################################
     ## set trace vars (this doesn't store the values, just updates it. ##
     ## Important in case user interrupts the burnin)                   ##
     #####################################################################

     s0_trace[[2]] <<- c(s0_trace[[2]], s0)
     s1_trace[[2]] <<- c(s1_trace[[2]], s1)
     tau_trace[[2]] <<- c(tau_trace[[2]], tau)
     alpha_r_trace[[2]] <<- rbind(alpha_r_trace[[2]], alpha_r)
     Yg_trace[[2]] <<- rbind(Yg_trace[[2]], Yg)
     variance_g_trace[[2]] <<- rbind(variance_g_trace[[2]], variance_g)

     inactive_means_trace[[2]] <<- c(inactive_means_trace[[2]], inactive_means)
     inactive_variances_trace[[2]] <<- c(inactive_variances_trace[[2]], inactive_variances)
     spike_probability_trace[[2]] <<- spike_probability
     active_means_trace[[2]] <<- rbind(active_means_trace[[2]], active_means)
     active_variances_trace[[2]] <<- rbind( active_variances_trace[[2]], active_variances)
     weight_active_trace[[2]] <<- weight_active
     weight_within_active_trace[[2]] <<- weight_within_active

     if(temperature == 1){

       all_allocation_burnin = allocation_active_inactive * allocation_within_active[[1]]

       allocation_trace <<- matrix(apply(component_matrix, 2, function(comp_matrix_col) 1 *
                                           (comp_matrix_col == all_allocation_burnin)), nrow = num_transcripts)

     }

     lnl_trace[[length(lnl_trace) + 1]] <<- .self$calculate_lnl(num_libraries)

     if(write_to_files & temperature == 1){

       .self$writeToOutputFiles(paste0(mcmc_prefixdir,"/", prefix), gen = i)

       if((i / sample_frequency) %% 4 == 0)
         .self$writeToYgVariancegOutputFiles(paste0(mcmc_prefixdir,"/", prefix), gen = i)

     }

     if(!is.null(target_ESS)  & length(lnl_trace) > 100){
       if(.self$calculate_lnl_ESS() > target_ESS){
         break
       }
     }

     cat("time: ", as.numeric(Sys.time()) - timestart, "\n")
     timestart = as.numeric(Sys.time())

   }

   i = i + 1

   gen <<- i

   ################################################
   ### tune the tuning parameters (tuningParam) ###
   ################################################
   if( gen %% 400 == 0 ) .self$tune_all(burnin_target_acceptance_rate)

  } #end mcmc

  lnl_trace <<- list(lnl_trace[[length(lnl_trace)]])

  ####################
  # set mirror moves #
  ####################
  mirror_moves <<- TRUE
  ntraceidx <- ceiling(0.75 * length(s0_trace[[2]])):length(s0_trace[[2]])

  s0_trace[[3]] <<- c(mean(s0_trace[[2]][ntraceidx]),
                     sd(s0_trace[[2]][ntraceidx]))
  s1_trace[[3]] <<- c(mean(log(-s1_trace[[2]])[ntraceidx]),
                      sd(log(-s1_trace[[2]])[ntraceidx]))
  tau_trace[[3]] <<- c(mean(log(tau_trace[[2]])[ntraceidx]),
                       sd(log(tau_trace[[2]])[ntraceidx]))
  alpha_r_trace[[3]] <<- rbind(colMeans(log(alpha_r_trace[[2]])[ntraceidx,]),
                               colSds(log(alpha_r_trace[[2]])[ntraceidx,]))
  Yg_trace[[3]] <<- cbind(colMeans(Yg_trace[[2]][ntraceidx,]), colSds(Yg_trace[[2]][ntraceidx,]))
  variance_g_trace[[3]] <<- cbind(colMeans(log(variance_g_trace[[2]][ntraceidx,])), colSds(log(variance_g_trace[[2]][ntraceidx,])))

  inactive_means_trace[[3]] <<- c(mean(log(-inactive_means_trace[[2]][ntraceidx])),
                                  sd(log(-inactive_means_trace[[2]][ntraceidx])))
  inactive_variances_trace[[3]] <<- c(mean(log(inactive_variances_trace[[2]][ntraceidx])),
                                            sd(log(inactive_variances_trace[[2]][ntraceidx])))

  active_means_dif_trace <- t(t(active_means_trace[[2]][ntraceidx,]) - threshold_a)
  active_means_trace[[3]] <<- rbind(colMeans(log(active_means_dif_trace)),
                                    colSds(log(active_means_dif_trace)))
  active_variances_trace[[3]] <<- rbind(colMeans(log(active_variances_trace[[2]][ntraceidx,])),
                                    colSds(log(active_variances_trace[[2]][ntraceidx,])))

  }


)
