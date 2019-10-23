zigzag$methods(

#' burnin
#' @name burnin
#' @description Run burnin and tune proposal size parameters for hieararchical bayesian mixture model
#' @usage burnin(sample_frequency = 10, ngen = 1000, burnin_target_acceptance_rate=0.23,
#' threads = 1, target_ESS = NULL, progress_plot = FALSE, write_to_files = TRUE,
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
  burnin = function(sample_frequency = 10, ngen = 1000, burnin_target_acceptance_rate=0.44,
                 threads = 1, target_ESS = NULL, progress_plot = FALSE, write_to_files = TRUE,
                 burninprefix = "output", append = FALSE){

  timestart = as.numeric(Sys.time())

  # initialize output files
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


  proposal_list <- list(.self$gibbsMixtureWeights, .self$gibbsAllocationActiveInactive, .self$gibbsAllocationWithinActive,
                       .self$mhInactiveMeans, .self$mhInactiveVariances,
                       .self$mhActiveMeansDif,
                       .self$mhActiveVariances,
                       .self$gibbsSpikeProb, .self$mhSpikeAllocation,
                       .self$mhYg, .self$mhSigma_g,
                       .self$mhTau, .self$mhSg, .self$mhS0Tau,
                       .self$mhP_x)

  plist_length <- length(proposal_list)

  ################
  # run the mcmc #
  ################

  # MCMC

  while(j <= ngen ){
   end = FALSE
   for(p in sample(seq(plist_length), plist_length, replace = TRUE, prob = proposal_probs)){
     proposal_list[[p]](tune = TRUE);

     if(is.na(sum(Yg)) || is.na(sum(sigma_g)) ||is.na(sum(tuningParam_yg)) || is.na(sum(inactive_spike_allocation)) ||
        is.nan(sum(inactive_spike_allocation)) || is.na(sum(allocation_active_inactive))||
        is.nan(sum(allocation_active_inactive)) || is.na(sum(sigma_g_trace)) || is.na(sum(Yg_trace))){

       print("an NA or NaN somewhere!")
       print(p)
       end = TRUE
       break

     }
   }
   if(end) break


   if(i %% sample_frequency == 0 & i > 2 * sample_frequency){

     xx_timestart = as.numeric(Sys.time())

     if(is.null(target_ESS)) j = i

     cat("#### ",i," ####  ", lnl_trace[[length(lnl_trace)]], "  prop. tscripts from inactive: ", round((sum(exp(Yg[inactive_idx])) - sum(exp(Yg[in_spike_idx])))/(10^6), 4),
         " inactive mean: ", round(inactive_means, 2), "  active_means: ", round(active_means, 2),"\n")

     ##########################
     ### Progress plot ########
     ##########################

     if(progress_plot){

       #plot variance trend
       plotdx = seq(-6,8,by=0.1)
       # plot(rowMeans(Xg),rowVars(Xg), xlim=c(-6,8), ylim=c(0,5),col=rgb(0,0,0,0.1))
       plot(Yg,sigma_g, xlim=c(-6,8), ylim = c(0, 0.5 * max(sigma_g)),col=rgb(0,0,0,0.1))
       polygon(c(plotdx, rev(plotdx)),
               c(qlnorm(0.025, s0 + s1*plotdx + tau, sqrt(tau)), rev(qlnorm(0.975, s0 + s1*plotdx +tau, sqrt(tau)))),
               col = rgb(1,0,0,0.2))
       lines(plotdx, exp(s0 + s1*plotdx),lwd = 2, col="red")

       #plot prob detection curves
       plot(seq(-10,10,by=0.1), .self$get_px(yy = seq(-10,10,by=0.1), gl = 1)[,1],type="l",ylim=c(0,1), ylab = "prob. detection", xlab = "log expression", col=rgb(0,0,0,0.1))
       sapply(1:num_libraries,function(x){lines(seq(-10,10,by=0.1), .self$get_px(yy = seq(-10,10,by=0.1), gl = 1)[,x])})

       #plot(mean Xg, Yg relationship)
       rawdat<-Xg; rawdat[rawdat == -Inf] <- -15
       plot(rowMeans(rawdat)[allocation_active_inactive==0],Yg[allocation_active_inactive==0],ylab = "Y", xlab = "log expression", col=rgb(0,0,1,0.2), xlim=c(-15,12),ylim=c(-8,12))
       points(rowMeans(rawdat)[allocation_active_inactive==1], Yg[allocation_active_inactive==1],col=rgb(1,0,0,0.2))
       abline(0,1,col="red")

       #plot Yg mixture distribution
       .self$plotHistoDensities(Yg)

     }

     #####################################################################
     ## set trace vars (this doesn't store the values, just updates it. ##
     ## Important in case user interrupts the burnin)                   ##
     #####################################################################

     s0_trace[[2]] <<- s0
     s1_trace[[2]] <<- s1
     alpha_r_trace[[2]] <<- alpha_r

     inactive_means_trace[[2]] <<- inactive_means
     inactive_variances_trace[[2]] <<- inactive_variances
     spike_probability_trace[[2]] <<- spike_probability
     active_means_trace[[2]] <<- active_means
     active_variances_trace[[2]] <<- active_variances
     weight_active_trace[[2]] <<- weight_active
     weight_within_active_trace[[2]] <<- weight_within_active

     if(temperature == 1){

       all_allocation_burnin = allocation_active_inactive * allocation_within_active[[1]]

       allocation_trace <<- matrix(apply(component_matrix, 2, function(comp_matrix_col) 1 *
                                           (comp_matrix_col == all_allocation_burnin)), nrow = num_transcripts)

     }

     lnl_trace[[length(lnl_trace) + 1]] <<- .self$calculate_lnl(num_libraries)

     if(write_to_files & temperature == 1) .self$writeToOutputFiles(paste0(mcmc_prefixdir,"/", prefix), gen = i)

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

   if( gen %% 400 == 0 ){

     tuningParam_alpha_r <<- sapply(1:num_libraries, function(xtrace){
       return(.self$x_tune(alpha_r_trace[[1]][[1]][xtrace,], tuningParam_alpha_r[xtrace], burnin_target_acceptance_rate, mintuningParam = 0.01, maxtuningParam = 2)
       )})
     tuningParam_s0 <<- .self$x_tune(s0_trace[[1]][[1]], tuningParam_s0, target_rate = burnin_target_acceptance_rate, mintuningParam = 0.0001, maxtuningParam = 5)
     tuningParam_s1 <<- .self$x_tune(s1_trace[[1]][[1]], tuningParam_s1, burnin_target_acceptance_rate, mintuningParam = 0.0001, maxtuningParam = 5)
     tuningParam_tau <<- .self$x_tune(tau_trace[[1]][[1]], tuningParam_tau, target_rate = burnin_target_acceptance_rate, mintuningParam = 0.0001, maxtuningParam = 5)
     tuningParam_s0tau <<- .self$x_tune(s0tau_trace[[1]][[1]], tuningParam_s0tau, burnin_target_acceptance_rate, mintuningParam = 0.0001, maxtuningParam = 5)
     tuningParam_sigma_g <<- .self$x_tune(sigma_g_trace, tuningParam_sigma_g, burnin_target_acceptance_rate, mintuningParam = 0.0001, maxtuningParam = 10)
     tuningParam_yg <<- .self$x_tune(Yg_trace, tuningParam_yg, burnin_target_acceptance_rate, mintuningParam = 0.0001, maxtuningParam = 10)

     tuningParam_multi_sigma <<- .self$x_tune(multi_sigma_trace[[1]][[1]], tuningParam_multi_sigma, burnin_target_acceptance_rate, mintuningParam = 0.001, maxtuningParam = 10)
     inactive_mean_tuningParam     <<- .self$x_tune(inactive_means_trace[[1]][[1]], inactive_mean_tuningParam, burnin_target_acceptance_rate, mintuningParam = 0.001, maxtuningParam = 10)
     inactive_variance_tuningParam <<- .self$x_tune(inactive_variances_trace[[1]][[1]], inactive_variance_tuningParam, burnin_target_acceptance_rate,
                                                    mintuningParam = 0.01, maxtuningParam = (inactive_variances_prior_log_max - inactive_variances_prior_log_min))

     active_mean_tuningParam       <<- sapply(seq(num_active_components), function(lib){
       return(.self$x_tune(active_means_trace[[1]][[1]][lib,], active_mean_tuningParam[lib], burnin_target_acceptance_rate, mintuningParam = 0.01, maxtuningParam = 10))})

     if(shared_active_variances){

       active_variance_tuningParam   <<- sapply(seq(num_active_components), function(lib){
         return(.self$x_tune(active_variances_trace[[1]][[1]][1,], active_variance_tuningParam[1], burnin_target_acceptance_rate,
                             mintuningParam = 0.01, maxtuningParam = (active_variances_prior_log_max - active_variances_prior_log_min)))})

     }else{

       active_variance_tuningParam   <<- sapply(seq(num_active_components), function(lib){
         return(.self$x_tune(active_variances_trace[[1]][[1]][lib,], active_variance_tuningParam[lib], burnin_target_acceptance_rate,
                             mintuningParam = 0.01, maxtuningParam = (active_variances_prior_log_max - active_variances_prior_log_min)))})

     }

     spike_probability_tuningParam <<- 1/.self$x_tune(spike_probability_trace[[1]][[1]], 1/spike_probability_tuningParam, burnin_target_acceptance_rate, maxtuningParam = 1/10, mintuningParam = 1/100000)
     mixture_weight_tuningParam <<- 1/.self$x_tune(mixture_weight_trace[[1]][[1]], 1/mixture_weight_tuningParam, burnin_target_acceptance_rate, maxtuningParam = 1/10, mintuningParam = 1/1000000)

   }

  } #end mcmc

  lnl_trace <<- list(lnl_trace[[length(lnl_trace)]])

  }


)
