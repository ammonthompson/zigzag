zigzag$methods(

  get_lib_bias_matrix = function(ib = inactive_bias, ab = active_bias, i_idx = inactive_idx, a_idx = active_idx){

    lb_matrix <- matrix(0, nrow = num_transcripts, ncol = num_libraries)
    for(lib in seq(num_libraries)){
      lb_matrix[i_idx,lib] <- ib[lib]
      lb_matrix[a_idx,lib] <- ab[lib]
     }

    return(lb_matrix)

  },

  trapz_integration = function(x, y){

    n = length(x)
    integral = 0.5 * sum((x[2:n] - x[1:(n-1)]) * (y[2:n] + y[1:(n-1)]))
    return(integral)

  },

  trapz_cumulative_integration = function(x, y){

    #gives cumulative integral at each x except x[1]
    n = length(x)
    integral = sapply(seq(n-1)+1, function(b) 0.5 * sum((x[2:b] - x[1:(b-1)]) * (y[2:b] + y[1:(b-1)])))
    return(integral)

  },

  get_model_density = function(expression_bins){

    #prior density of Yg out of spike
    ds =  ((1 - weight_active) * (1 - spike_probability) * dnorm(expression_bins, inactive_means, sqrt(inactive_variances)) +
             weight_active * rowSums( sapply(seq(length(weight_within_active)), function(k){
               weight_within_active[k] * dnorm(expression_bins, active_means[k], sqrt(active_variances[k]))
             })))/((1 - weight_active) * (1 - spike_probability) + weight_active)

    return(ds)
  },

  get_model_cumulative = function(expression_bins){
    #returns approximate cumulative distribution at bin values for model excluding bin[1]

    ds = get_model_density(expression_bins)

    cum_ds = trapz_cumulative_integration(expression_bins, ds)

    return(cum_ds)

  },

  r_dirichlet = function(alpha){

    mygamma <- rgamma(length(alpha), shape = alpha, rate = 1)

    return(mygamma/sum(mygamma))

  },

  d_dirichlet = function(x, alpha, log = FALSE) {

    if (is.matrix(x)) return( apply(x, 1, d_dirichlet) )

    if( any(x < 0)  ) return(0L)

    tol <- .Machine$double.eps

    if( abs(sum(x) - 1) > tol ) return(0L)

    stopifnot(length(x) == length(alpha))

    log_f <- sum((alpha - 1)*log(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha))

    if(log) return(log_f)

    exp(log_f)

  },

  calculate_lnl = function(num_libs){

    if(num_libs > 1){

      yg_lnl <- sum(.self$computeXgLikelihood(Xg, Yg, variance_g, p_x))

    }else{

      yg_lnl <- sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx])) +
              sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], active_means, active_variances))

    }

    return(yg_lnl)

  },

  calculate_lnl_ESS = function(){

    lnltrace = unlist(lnl_trace)[-1]

    ess <- as.numeric(coda::effectiveSize(as.mcmc(lnltrace)))

    if(ess < length(lnltrace)){

      return(ess)

    }else{

      return(length(lnltrace)/(1 + 2 * sum(cor(lnltrace[-1], lnltrace[-length(lnltrace)]))))

    }


  },

  calculate_active_means = function(am = active_means_dif, mt = multi_thresh_a){

    if(mt){

      return(am + threshold_a)

    }else{

      return(cumsum(am) + threshold_a)

    }

  },

  setActiveInactive_idx = function(){
    active_idx <<- which(allocation_active_inactive == 1)
    inactive_idx <<- which(allocation_active_inactive == 0)
  },

  setInSpike_idx = function(){

    in_spike_idx <<- which(inactive_spike_allocation == 1)

    out_spike_idx <<- which(inactive_spike_allocation == 0)

  },

  set_sigmaX_pX = function(){

    Sg <<- exp(s0 + s1 * Yg)

    # p_x <<- matrix(sapply(1:num_libraries, function(lib){return(1 - exp(-alpha_r[lib] * gene_lengths * exp(Yg)))}),
    #                nrow = num_transcripts)
    # p_x <<- 1 - exp(-outer((gene_lengths * exp(Yg)), alpha_r))
    p_x <<- .self$get_px()

  },

  get_sigmax = function(ss0 = s0, ss1 = s1, yy = Yg){

    return(exp(ss0 + ss1 * yy))

  },

  get_px = function(falpha_r = alpha_r, yy = Yg, gl = gene_lengths, recover_x = F, bias_m = lib_bias_matrix){

    if(recover_x) recover()
    # return(1 - exp(-outer((gl * exp(yy)), falpha_r)))
    return(1 - exp(-(gl * exp(yy + bias_m)) * matrix(falpha_r, nrow = length(yy), ncol = length(falpha_r), byrow = T)))

  },

  x_tune = function(param_trace_vector, tuningParam, target_rate = 0.44, mintuningParam = 0, maxtuningParam = Inf){

    if(is.vector(param_trace_vector)){

      pjump = sum(param_trace_vector)/100

    }else{

      pjump = rowSums(param_trace_vector)/100

    }

    new_tuningParam = tuningParam * (tan(pi * pjump/2)/tan(pi * target_rate/2))

    new_tuningParam[which(new_tuningParam < mintuningParam)] = mintuningParam

    new_tuningParam[which(new_tuningParam > maxtuningParam)] = maxtuningParam

    return(new_tuningParam)

  },

  tune_all = function(burnin_target_acceptance_rate){

    tuningParam_alpha_r <<- sapply(1:num_libraries, function(xtrace){
      return(.self$x_tune(alpha_r_trace[[1]][[1]][xtrace,], tuningParam_alpha_r[xtrace],
                          burnin_target_acceptance_rate, mintuningParam = 0.01, maxtuningParam = 2)
      )})

    tuningParam_s0 <<- .self$x_tune(s0_trace[[1]][[1]], tuningParam_s0,
                                    target_rate = burnin_target_acceptance_rate,
                                    mintuningParam = 0.0001, maxtuningParam = 5)
    tuningParam_s1 <<- .self$x_tune(s1_trace[[1]][[1]], tuningParam_s1,
                                    burnin_target_acceptance_rate,
                                    mintuningParam = 0.0001, maxtuningParam = 5)
    tuningParam_tau <<- .self$x_tune(tau_trace[[1]][[1]], tuningParam_tau,
                                     target_rate = burnin_target_acceptance_rate,
                                     mintuningParam = 0.0001,  maxtuningParam = 5)
    tuningParam_s0tau <<- .self$x_tune(s0tau_trace[[1]][[1]], tuningParam_s0tau,
                                       burnin_target_acceptance_rate,
                                       mintuningParam = 0.0001, maxtuningParam = 5)
    tuningParam_variance_g <<- .self$x_tune(variance_g_trace[[1]], tuningParam_variance_g,
                                         burnin_target_acceptance_rate,
                                         mintuningParam = 0.0001, maxtuningParam = 100)
    tuningParam_yg <<- .self$x_tune(Yg_trace[[1]], tuningParam_yg, burnin_target_acceptance_rate,
                                    mintuningParam = 0.0001, maxtuningParam = 100)

    tuningParam_multi_sigma <<- .self$x_tune(multi_sigma_trace[[1]][[1]], tuningParam_multi_sigma,
                                             burnin_target_acceptance_rate,
                                             mintuningParam = 0.001, maxtuningParam = 10)
    inactive_mean_tuningParam     <<- .self$x_tune(inactive_means_trace[[1]][[1]], inactive_mean_tuningParam,
                                                   burnin_target_acceptance_rate,
                                                   mintuningParam = 0.001, maxtuningParam = 10)
    inactive_variance_tuningParam <<- .self$x_tune(inactive_variances_trace[[1]][[1]], inactive_variance_tuningParam,
                                                   burnin_target_acceptance_rate,
                                                   mintuningParam = 0.01, maxtuningParam = (inactive_variances_prior_log_max - inactive_variances_prior_log_min))

    active_mean_tuningParam       <<- sapply(seq(num_active_components), function(lib){
      return(.self$x_tune(active_means_trace[[1]][[1]][lib,], active_mean_tuningParam[lib],
                          burnin_target_acceptance_rate,
                          mintuningParam = 0.01, maxtuningParam = 10))})

    if(shared_active_variances){

      active_variance_tuningParam   <<- sapply(seq(num_active_components), function(lib){
        return(.self$x_tune(active_variances_trace[[1]][[1]][1,], active_variance_tuningParam[1], burnin_target_acceptance_rate,
                            mintuningParam = 0.01, maxtuningParam = (active_variances_prior_log_max - active_variances_prior_log_min)))})

    }else{

      active_variance_tuningParam   <<- sapply(seq(num_active_components), function(lib){
        return(.self$x_tune(active_variances_trace[[1]][[1]][lib,], active_variance_tuningParam[lib], burnin_target_acceptance_rate,
                            mintuningParam = 0.01, maxtuningParam = (active_variances_prior_log_max - active_variances_prior_log_min)))})

    }

    spike_probability_tuningParam <<- 1/.self$x_tune(spike_probability_trace[[1]][[1]], 1/spike_probability_tuningParam, burnin_target_acceptance_rate, maxtuningParam = 1/10, mintuningParam = 1/10000000)
    mixture_weight_tuningParam <<- 1/.self$x_tune(mixture_weight_trace[[1]][[1]], 1/mixture_weight_tuningParam, burnin_target_acceptance_rate, maxtuningParam = 1/10, mintuningParam = 1/100000000)

    tuningParam_inactive_bias <<- 1/.self$x_tune(inactive_bias_trace[[1]][[1]], 1/tuningParam_inactive_bias, target_rate = burnin_target_acceptance_rate/2,
                                                 maxtuningParam = 1/10, mintuningParam = 1/1000000000)
    tuningParam_active_bias <<- 1/.self$x_tune(active_bias_trace[[1]][[1]], 1/tuningParam_active_bias, target_rate = burnin_target_acceptance_rate/2,
                                               maxtuningParam = 1/10, mintuningParam = 1/1000000000)

#
#     # print(tuningParam_variance_g[1:10])
#     write.table(matrix(tuningParam_variance_g, nrow = 1), file = "check_variance_g_tuningparams.txt", sep = "\t", quote = F, append = T, row.names = F, col.names = F )

  },

  get_logPosterior = function(t = temperature){

    if(num_libraries > 1){

      lnLxy <- sum(.self$computeXgLikelihood(Xg, Yg, Sg, p_x)) +
        sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx])) +
        sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], active_means, active_variances))

      lnPrior <- .self$computeSpikePriorProbability(spike_probability) + .self$computeAlphaRPriorProbability(alpha_r) +
        .self$computeActiveWeightPriorProbability(weight_active) + .self$computeWithinActiveWeightPriorProbability(weight_within_active) +
        .self$computeActiveMeansDifPriorProbability(active_means_dif) + .self$computeActiveVariancesPriorProbability(active_variances) +
        .self$computeInactiveMeansPriorProbability(inactive_means) + .self$computeInactiveVariancesPriorProbability(inactive_variances) +
        .self$computeVarianceGPriorProbability(variance_g, Sg) + .self$computeTauPriorProbability(tau) + .self$computeS0PriorProbability(s0) + .self$computeS1PriorProbability(s1)

      tlogpost <- t * sum(lnLxy + lnPrior)

    }else{

      lnLxy <- sum(.self$computeInactiveLikelihood(Yg[inactive_idx], spike_probability, inactive_means, inactive_variances, inactive_spike_allocation[inactive_idx])) +
        sum(.self$computeActiveLikelihood(Yg[active_idx], allocation_within_active[[1]][active_idx], active_means, active_variances))

      lnPrior <- .self$computeSpikePriorProbability(spike_probability) +
        .self$computeActiveWeightPriorProbability(weight_active) + .self$computeWithinActiveWeightPriorProbability(weight_within_active) +
        .self$computeActiveMeansDifPriorProbability(active_means_dif) + .self$computeActiveVariancesPriorProbability(active_variances) +
        .self$computeInactiveMeansPriorProbability(inactive_means) + .self$computeInactiveVariancesPriorProbability(inactive_variances)

      tlogpost <- t * (lnLxy + lnPrior)

    }

    return(tlogpost)

  },

  get_fields = function(){

    return(field_list_values)

  },

  set_fields = function(f_list){  # So far doesn't work

    for(i in length(f_list)){

      assign(field_list_names[[i]], f_list[[i]], envir=.self)

      # field(field_list_names[[i]], f_list[[i]])

    }

  },

  get_tuningParam_fields = function(){

    tuningParam_list <- list(inactive_mean_tuningParam,
                             inactive_variance_tuningParam,
                             spike_probability_tuningParam,
                             active_mean_tuningParam,
                             active_variance_tuningParam,
                             mixture_weight_tuningParam,
                             tuningParam_s0,
                             tuningParam_s1,
                             tuningParam_tau,
                             tuningParam_s0tau,
                             tuningParam_alpha_r,
                             tuningParam_yg,
                             tuningParam_variance_g,
                             tuningParam_multi_sigma,
                             tuningParam_sigma_mu)


    return(tuningParam_list)

  },

  set_tuningParam_fields = function(hlist){


    inactive_mean_tuningParam <<- hlist[[1]]
    inactive_variance_tuningParam <<- hlist[[2]]
    spike_probability_tuningParam <<- hlist[[3]]
    active_mean_tuningParam <<- hlist[[4]]
    active_variance_tuningParam <<- hlist[[5]]
    mixture_weight_tuningParam <<- hlist[[6]]
    tuningParam_s0 <<- hlist[[7]]
    tuningParam_s1 <<- hlist[[8]]
    tuningParam_tau <<- hlist[[9]]
    tuningParam_s0tau <<- hlist[[10]]
    tuningParam_alpha_r <<- hlist[[11]]
    tuningParam_yg <<- hlist[[12]]
    tuningParam_variance_g <<- hlist[[13]]
    tuningParam_multi_sigma <<- hlist[[14]]
    tuningParam_sigma_mu <<- hlist[[15]]


  },

  set_temperature = function(tt){

    temperature <<- tt

  },

  ###################################
  #### TESTING or DEPRICATED ########
  ###################################

  xxget_px = function(falpha_r = alpha_r, yy = Yg, sg = variance_g, xx = Xg, gl = gene_lengths, recover_x = FALSE){

    # doesn't work. Not sure why, but when mhP_x is run, XgLikelihood seems to increment more and more positive for a lot of genes.

    sd_variance_g = sqrt(sg)

    if(recover_x == TRUE) recover()

    if(is.matrix(xx)){

      npdetect = sapply(1:length(falpha_r), function(lib){

        yyzero = which(xx[,lib] == -Inf)

        yynotzero = which(xx[,lib] > -Inf)

        lib_pdetect = NULL

        #P(x is not detected)
        lib_pdetect[yyzero] = 1 - as.numeric(sapply(yyzero, function(gidx){
          integrate(function(x){exp(-falpha_r[lib] * gl[gidx] * exp(x)) * dnorm(x, yy[gidx], sd_variance_g[gidx])},
                    lower = yy[gidx] - 8*sd_variance_g[gidx], upper = yy[gidx] + 8*sd_variance_g[gidx])[1]
        }))


        #p(x is detected)
        lib_pdetect[yynotzero] = 1 - exp(-falpha_r[lib] * gl[yynotzero] * exp(xx[yynotzero,lib]))

        return(lib_pdetect)

      })

      return(npdetect)


    }else{

      yyzero = which(xx == -Inf)

      yynotzero = which(xx > -Inf)

      npdetect = NULL

      #P(x is not detected)
      npdetect[yyzero] = 1 - as.numeric(sapply(yyzero, function(gidx){
        integrate(function(x){exp(-falpha_r * gl[gidx] * exp(x)) * dnorm(x, yy[gidx], sd_variance_g[gidx])},
                  lower = yy[gidx] - 8*sd_variance_g[gidx], upper = yy[gidx] + 8*sd_variance_g[gidx])[1]
      }))


      #P(x is not detected)
      npdetect[yynotzero] = 1 - exp(-falpha_r * gl[yynotzero] * exp(xx[yynotzero]))

      return(npdetect)

    }


  },


  test_tune = function(param, tuningParam, active_param = FALSE, target_rate, maxtuningParam, inverse = FALSE){

    sfpower = 1

    if(inverse) sfpower = -1

    if(active_param){

      param_rate=lapply(param[[1]],function(x){return(sapply(1:num_active_components,function(y){return(sum(x[y,])/100)}))}) # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

      scale_factor=lapply(param_rate,function(x){
        return(sapply(1:num_active_components, function(y){
          return(exp(sfpower * 0.05 * (x[y] - target_rate)/target_rate))
        }))
      })

      newtuningParam_list <- sapply(1:num_active_components, function(y){
        newtuningParam = scale_factor[[1]][y]*tuningParam[y]

        if(newtuningParam > maxtuningParam){
          return(maxtuningParam)
        }else{
          return(newtuningParam)
        }

      })

    }else{


      if(is.list(param)){

        param_rate=lapply(param[[1]],function(x){return(sum(x)/100)}) # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

        scale_factor=lapply(param_rate,function(x){
          return(exp(sfpower * 0.05 * (x - target_rate)/target_rate))
        })


        newtuningParam <- sapply(1:length(tuningParam), function(x){
          newtuningParam = scale_factor[[x]] * tuningParam[x]

          if(newtuningParam > maxtuningParam){
            return(maxtuningParam)
          }else{
            return(newtuningParam)
          }

        })

      }else{

        param_rate=sum(param)/100 # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

        scale_factor=lapply(param_rate,function(x){

          return(exp(sfpower * 0.05 * (x - target_rate)/target_rate))

        })

        newtuningParam_list <- sapply(1:length(tuningParam), function(x){
          newtuningParam = scale_factor[[x]]*tuningParam[x]

          if(newtuningParam > maxtuningParam){
            return(maxtuningParam)
          }else{
            return(newtuningParam)
          }
        })

      }
      recover()
    }
    return(newtuningParam_list)

  },

  tune_logit = function(param, tuningParam, active_param = FALSE, target_rate, maxtuningParam, inverse = FALSE){

    sfpower = 1

    if(inverse) sfpower = -1


    param_rate=lapply(param[[1]],function(x){return(sum(x)/100)}) # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}


    scale_factor=lapply(param_rate,function(x){

      return(exp(sfpower * 0.05 * (x - target_rate)/target_rate))

    })

    newtuningParam_list <- sapply(1:length(tuningParam), function(x){
      newtuningParam = scale_factor[[x]]*tuningParam[x]

      if(newtuningParam > maxtuningParam){
        return(maxtuningParam)
      }else{
        return(newtuningParam)
      }
    })


  },

  y_tune = function(param_trace_vector, tuningParam, target_rate = 0.44, mintuningParam = 0, maxtuningParam = Inf){

    if(is.vector(param_trace_vector)){

      pjump = sum(param_trace_vector)/100

    }else{

      pjump = rowSums(param_trace_vector)/100

    }

    # double p = this->targetAcceptanceRate;
    # if ( rate > p )
    # {
    #   alpha_r *= (1.0 + ((rate-p)/(1.0 - p)) );
    # }
    # else
    # {
    #   alpha_r /= (2.0 - rate/p);
    # }

    if(pjump > target_rate){

      new_tuningParam = tuningParam * (1 + ((pjump - target_rate)/(1 - target_rate)))

    }else{

      new_tuningParam = tuningParam * target_rate/(2 - pjump)
    }

    new_tuningParam[which(new_tuningParam < mintuningParam)] = mintuningParam

    new_tuningParam[which(new_tuningParam > maxtuningParam)] = maxtuningParam

    return(new_tuningParam)

  },

  tune = function(param, tuningParam, active_param = FALSE, target_rate, maxtuningParam, inverse = FALSE){

    change_scale = 0.1

    sfpower = 1

    if(inverse) sfpower = -1

    if(active_param){

      param_rate=lapply(param[[1]],function(x){return(sapply(1:num_active_components,function(y){return(sum(x[y,])/100)}))}) # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

      scale_factor=lapply(param_rate,function(x){

        return(sapply(1:num_active_components, function(y){
          return(exp(sfpower * change_scale * (x[y] - target_rate)/target_rate))
        }))

      })

      newtuningParam_list <- sapply(1:num_active_components, function(y){
        newtuningParam = scale_factor[[1]][y]*tuningParam[y]

        if(newtuningParam > maxtuningParam){
          return(maxtuningParam)
        }else{
          return(newtuningParam)
        }

      })

    }else{

      param_rate=lapply(param[[1]],function(x){return(sum(x)/100)}) # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

      if(is.list(param)){

        param_rate=lapply(param[[1]],function(x){return(sum(x)/100)}) # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

        scale_factor=lapply(param_rate,function(x){
          return(exp(sfpower * change_scale * (x - target_rate)/target_rate))
        })


        newtuningParam_list <- sapply(1:length(tuningParam), function(x){
          newtuningParam = scale_factor[[x]]*tuningParam[x]

          if(newtuningParam > maxtuningParam){
            return(maxtuningParam)
          }else{
            return(newtuningParam)
          }

        })

      }else{

        param_rate=sum(param)/100 # [[1]] is a record of the previous 100 proposals and if they were accepted or not {0, 1}

        scale_factor=lapply(param_rate,function(x){

          return(exp(sfpower * change_scale * (x - target_rate)/target_rate))

        })

        newtuningParam_list <- sapply(1:length(tuningParam), function(x){
          newtuningParam = scale_factor[[x]]*tuningParam[x]

          if(newtuningParam > maxtuningParam){
            return(maxtuningParam)
          }else{
            return(newtuningParam)
          }
        })

      }

    }
    return(newtuningParam_list)

  } #depricated

)
