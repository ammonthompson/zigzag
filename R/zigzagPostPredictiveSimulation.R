zigzag$methods(

  posteriorPredictiveSimulation = function(num_sims = 1000, plotlib = 1, prefix = "out", recover_x = FALSE){

    pp_prefix = paste0(output_directory, "/", prefix, "_mcmc_output/", prefix)

    #post pred statistics pars
    binwidth = 0.1
    model_range = c(qnorm(0.000001, inactive_means, sqrt(inactive_variances)),
                    qnorm(0.999999, active_means[num_active_components], sqrt(active_variances[num_active_components])))
    Yg_range = c(min(Yg[out_spike_idx]), max(Yg[out_spike_idx]))
    xg_greater_minusInf_idx = which(Xg > -Inf, arr.ind = T)
    Xg_range = c(min(Xg[xg_greater_minusInf_idx]), max(Xg[xg_greater_minusInf_idx]))
    all_ranges = c(model_range, Yg_range, Xg_range)

    bins = seq(min(all_ranges) - 2, max(all_ranges) + 2, by = binwidth)
    yg_lessThan_bin = sapply(bins, function(y) (Yg < y) * 1 )


    #### Simulate data from posterior and compute realized discrepancy statstics

    sim_xg = sapply(1:num_libraries, function(x){return(rnorm(length(Yg), Yg, sqrt(variance_g)))})

    # p detect
    sim_xg = sapply(1:num_libraries, function(lib){

      temp_sim_xg <- sim_xg[,lib]
      not_detected <- (runif(length(Yg)) < (1 - p_x[,lib]))
      temp_sim_xg[not_detected] <- -Inf

      return(temp_sim_xg)

    })

    ### simulate spike p( spike | inactive )
    sim_xg[in_spike_idx,] <- rep(-Inf, num_libraries)

    sim_xg <- as.data.frame(sim_xg)

    # lower level discrepencies, zeros and non-zero Xg distros
    W_L <- .self$get_lowerLevel_wasserstein(Xg, bins, yg_lessThan_bin) -
           .self$get_lowerLevel_wasserstein(sim_xg, bins, yg_lessThan_bin)

    B_rumsfeld <- .self$get_rumsfeld(Xg) - .self$get_rumsfeld(sim_xg)


    if(recover_x) recover()

    # scaled library-specific wasserstein discrepancies. Xg scaled by common mean
    # so the discrepancy ignores library-specific scaling factor differences
    # which are captured by the gene-wise variance model.
    # if the model estimate library-specific biases for inactive and active components
    # then do not rescale real and simulated data
    if(library_bias){
      scaled_xg = Xg
    }else{
      lib_means_xg <- apply(Xg, 2, function(sx) mean(sx[sx > -Inf]))
      grand_mean_xg <- mean(lib_means_xg)
      scaled_xg <- sapply(seq(num_libraries), function(sx) Xg[,sx] - lib_means_xg[sx] + grand_mean_xg)
    }
    if(library_bias){
      scaled_sim_xg <- sim_xg
    }else{
      lib_means_sim_xg <- apply(sim_xg, 2, function(sx) mean(sx[sx > -Inf]))
      grand_mean_sim_xg <- mean(lib_means_sim_xg)
      scaled_sim_xg <- sapply(seq(num_libraries), function(sx) sim_xg[,sx] - lib_means_sim_xg[sx] + grand_mean_sim_xg)
    }

    W_L_r <- sapply(seq(ncol(scaled_xg)), function(sxg){

      .self$get_lowerLevel_wasserstein(scaled_xg[,sxg, drop = F], bins, yg_lessThan_bin, p_x[,sxg, drop = F]) -
        .self$get_lowerLevel_wasserstein(scaled_sim_xg[,sxg, drop = F], bins, yg_lessThan_bin, p_x[,sxg, drop = F])

    })


    ### upper level discrepency simulating Yg
    #simulate Yg
    sim_yg <- NULL
    sim_yg[inactive_idx] <- rnorm(length(inactive_idx), inactive_means, sqrt(inactive_variances))
    sim_yg[active_idx] <- rnorm(length(active_idx), active_means[allocation_within_active[[1]][active_idx]],
                               sqrt(active_variances[allocation_within_active[[1]][active_idx]]))


    W_U <- .self$get_upperLevel_wasserstein(Yg, bins) -
           .self$get_upperLevel_wasserstein(sim_yg, bins)

    cat("W_L: ",W_L, "  W_U: ", W_U, '\n')

    write.table(matrix(c(gen, W_L, W_U, B_rumsfeld[1]), nrow = 1),
                file = paste0(pp_prefix, ".post_predictive_output.log"),
                append=T, sep="\t", row.names=F, col.names=F)

    write.table(matrix(c(gen, W_L_r, B_rumsfeld[-1]), nrow = 1),
                file = paste0(pp_prefix, ".library_specific_post_predictive_output.log"),
                append=T, sep="\t", row.names=F, col.names=F)

    return(list(sim_xg, sim_yg))

  },

  get_cumulative = function(xdat, xbins){ #xdat is a vector (one library)

    zeroIdx = which(xdat == -Inf)
    num_zeros = length(zeroIdx)
    p_zeros = num_zeros/num_transcripts

    xdat = xdat[-zeroIdx]

    F = c()
    for( x in xbins ) F = c( F, matrixStats::count( xdat <= x )/(num_transcripts - num_zeros))

    return( c(p_zeros, F) )

  },

  get_rumsfeld = function(xdat){

    model_Br = sapply(seq(num_libraries), function(r) (1 - spike_probability * (1 - weight_active)) * mean(1 - p_x[out_spike_idx, r]) +
                        spike_probability * (1 - weight_active))

    xdat_Br = sapply(seq(num_libraries), function(r) sum(xdat[,r] == -Inf)/num_transcripts)

    delta_Br <- abs(xdat_Br - model_Br)

    return(c(mean(delta_Br), delta_Br))

  },

  get_lowerLevel_wasserstein = function(xdat, xbins, yg_lessThan_bin_infun, xp_x = p_x){

    xdat_num_libs <- ncol(xdat)

    model_cum = sapply(seq(xdat_num_libs), function(r){
      colMeans(xp_x[out_spike_idx, r] * yg_lessThan_bin_infun[out_spike_idx, ])
    })


    norm_factor = matrixStats::colMaxs(model_cum)

    model_cum = t(t(model_cum)/norm_factor)

    xdat_cum = sapply(seq(xdat_num_libs), function(r){
      sapply(xbins, function(y) mean(xdat[xdat[,r] > -Inf, r] < y))
    })

    #average number of genes detected < x averaged over all libraries and the wasserstein distance.
    mean_xdat_cum = rowMeans(xdat_cum)
    mean_model_cum = rowMeans(model_cum)
    Wass_means = trapz_integration(xbins, abs(mean_xdat_cum - mean_model_cum))
    return(c(Wass_means))

  },

  get_upperLevel_wasserstein = function(ydat, xbins){

    model_cum = get_model_cumulative(xbins)
    ydat_cum = sapply(xbins[-1], function(x) mean(ydat[out_spike_idx] < x))

    Wass_means = trapz_integration(xbins[-1], abs(ydat_cum - model_cum))

    return(Wass_means)

  },


  #postpred plot functions

  initialize_postPredPlots = function(multi_plot_pars_f, mcmc_prefixdir_f, prefix_f,
                                      multiL1PlotDevice, upperPlotDevice){},

  update_postPredPlots = function(post_pred_instance_f, multi_plot_pars_f,
                                  multiL1PlotDevice, upperPlotDevice, d_posXg_rowMeans_f){

      cat("plotting post pred...\n")

      ### Multiple libs plots ###
      num_libs_plot_postpred = num_libraries
      if(num_libs_plot_postpred > 25) num_libs_plot_postpred = 25

      dev.set(multiL1PlotDevice)

      lib_means_sim_xg = apply(post_pred_instance_f[[1]], 2, function(sx) mean(sx[sx > -Inf]))
      grand_mean_sim_xg = mean(lib_means_sim_xg)
      scaled_sim_xg = sapply(seq(num_libraries), function(sx) post_pred_instance_f[[1]][,sx] - lib_means_sim_xg[sx] + grand_mean_sim_xg)

      x_from <- min(Xg[which(Xg > inf_tol, arr.ind = T)]) * 0.99
      x_to <- max(Xg[which(Xg > inf_tol, arr.ind = T)]) * 1.01
      ppred_bw <- density(Xg[,1], from = x_from, to = x_to)$bw
      for(lib_plot in seq(num_libs_plot_postpred)){
        par("usr" = multi_plot_pars_f[[lib_plot]]$usr)
        par("plt" = multi_plot_pars_f[[lib_plot]]$plt)
        par("mfg" = multi_plot_pars_f[[lib_plot]]$mfg)

        d_inactive = density(scaled_sim_xg[inactive_idx,lib_plot], bw = ppred_bw, from = x_from, to = x_to)
        d_active = density(scaled_sim_xg[active_idx,lib_plot], bw = ppred_bw, from = x_from, to = x_to)

        # lines(density(scaled_sim_xg[,lib_plot]), col = "chartreuse3", lwd = 0.25)
        lines(d_inactive$x, d_inactive$y * (1 - weight_active) + d_active$y * weight_active, col = "chartreuse3", lwd = 0.25)

        lines(d_inactive$x, d_inactive$y * (1 - weight_active),col = "dodgerblue", lwd = 0.5)
        lines(d_active$x, d_active$y * weight_active, col = "tomato", lwd = 0.5)

        # p zero
        pzero = sum(post_pred_instance_f[[1]][,lib_plot] == inf_tol)/num_transcripts
        pzero_xpos = par("usr")[1] + 0.05 * diff(par("usr")[1:2])
        pzero_ypos = par("usr")[4]
        points(pzero_xpos, pzero * pzero_ypos, col = "dodgerblue")
      }


      ### Upper level plot ###
      dev.set(upperPlotDevice)

      par(mfg = c(1,1))

      # ybins = seq(min(c(Yg[out_spike_idx] , post_pred_instance_f[[2]][out_spike_idx])) -1 ,
      #             max(c(Yg[out_spike_idx] , post_pred_instance_f[[2]][out_spike_idx]))+1, by = 0.2)
      #
      # dYg = hist(Yg[out_spike_idx], plot = F, breaks = ybins)
      # dsimYg = hist(post_pred_instance_f[[2]][out_spike_idx], plot = F, breaks = ybins)

      lines(density(post_pred_instance_f[[2]][out_spike_idx], from = x_from, to = x_to, bw = ppred_bw), col = "chartreuse3", lwd =0.25)
      lines(density(Yg[out_spike_idx], from = x_from, to = x_to, bw = ppred_bw), col = "orange", lwd =0.25)

      par(mfg = c(2,1))
      dY_inactive = density(Yg[intersect(out_spike_idx, inactive_idx)], from = x_from, to = x_to, bw = ppred_bw)
      dY_active = density(Yg[active_idx], from = x_from, to = x_to, bw = ppred_bw)
      lines(density(Yg[out_spike_idx], from = x_from, to = x_to, bw = ppred_bw), col = "orange", lwd =0.25)
      lines(dY_inactive$x, dY_inactive$y * (1 - weight_active * num_transcripts/length(out_spike_idx)), col = "dodgerblue", lwd =0.25)
      lines(dY_active$x, dY_active$y * weight_active * num_transcripts/length(out_spike_idx), col = "tomato", lwd =0.25)

  },

  finalize_postPredPlots = function(multi_plot_pars_f,
                                    multiL1PlotDevice, upperPlotDevice, d_posXg_rowMeans_f, scaled_xg_f){

    num_libs_plot_postpred = num_libraries
    if(num_libs_plot_postpred > 25) num_libs_plot_postpred = 25

    #Level 1
    dev.set(multiL1PlotDevice)
    for(i in seq(num_libs_plot_postpred)){
      par("usr" = multi_plot_pars_f[[i]]$usr)
      par("plt" = multi_plot_pars_f[[i]]$plt)
      par("mfg" = multi_plot_pars_f[[i]]$mfg)
      lines(density(scaled_xg_f[,i]), lwd = 1)

      # p zero
      pzero = sum(Xg[,i] == inf_tol)/num_transcripts
      pzero_xpos = par("usr")[1] + 0.07 * diff(par("usr")[1:2])
      pzero_ypos = par("usr")[4]
      points(pzero_xpos, pzero * pzero_ypos, col = "black", pch = 4)

    }

    #Upper level (2)
    dev.set(upperPlotDevice)
    abline(h= -0.6 * max(d_posXg_rowMeans_f$y), lwd = 1.5)


    dev.off(multiL1PlotDevice)
    dev.off(upperPlotDevice)
  }

)
