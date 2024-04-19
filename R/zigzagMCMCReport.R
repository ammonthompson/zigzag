zigzag$methods(

  #Utilities for assessing and graphing mcmc output files

  generate_MCMC_reports = function(mcmc_dir, mcmc_file_prefix){

    mcmcReport_output_directory <- "mcmc_report"

    dir.create(paste0(mcmc_dir, "/", mcmcReport_output_directory))

    prfx <- paste0(mcmc_dir, "/", mcmc_file_prefix)

    cat("Analyzing MCMC files. Writing to files located in: ", mcmcReport_output_directory, "...\n")

    #load mcmc log files
    model_log <- read.table(paste0(prfx, "_model_parameters.log"), header = T, row.names = 1)
    yg_log <- read.table(paste0(prfx, "_yg_candidate_genes.log"), header = T, row.names = 1)
    varg_log <- read.table(paste0(prfx, "_varianceg_candidate_genes.log"), header = T, row.names = 1)

    mcmc_report_prefix = paste0(mcmc_dir, "/", mcmcReport_output_directory, "/", mcmc_file_prefix)

    alpha_r_idx <- grep(pattern = "library", colnames(model_log), value = FALSE)

    #create reports and plots
    .self$makePosteriorPriorPlots(model_log, paste0(mcmc_report_prefix, "_posterior_distributions.pdf"))
    .self$makeEssReport(model_log, yg_log, varg_log, mcmc_report_prefix)
    .self$makeMCMCplots(model_log[,-alpha_r_idx], paste0(mcmc_report_prefix, "_model_params_trace.pdf"))
    .self$makeMCMCplots(model_log[,alpha_r_idx], paste0(mcmc_report_prefix, "_library_params_trace.pdf"))

  },

  makeEssReport = function(param_mcmc_output, yg_mcmc_output, varg_mcmc_output, essPrefix){

    paramESS <- round(coda::effectiveSize(param_mcmc_output), digits = 2)

    ygESS <- round(coda::effectiveSize(yg_mcmc_output), digits = 2)

    vargESS <- round(coda::effectiveSize(varg_mcmc_output), digits = 2)

    warn <- c(paramESS[which(paramESS < 200)], ygESS[which(ygESS < 100)], vargESS[vargESS < 100])

    write.table(matrix(c("param", "ESS"), nrow = 1), file = paste0(essPrefix, "_model_ess.tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(paramESS, file = paste0(essPrefix, "_model_ess.tsv"),
                row.names = names(paramESS), col.names = F, append = TRUE,
                sep = "\t", quote = F)

    write.table(matrix(c("gene", "yg_ESS", "varg_ESS"), nrow = 1), file = paste0(essPrefix, "_Yg_varg_ess.tsv"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(cbind(ygESS, vargESS), file = paste0(essPrefix, "_Yg_varg_ess.tsv"),
                row.names = names(ygESS), col.names = F, quote = F, sep = "\t", append = TRUE)

    write.table(matrix(c("param_or_gene", "ESS"), nrow = 1), file = paste0(essPrefix, "_WARNINGS_ess.txt"),
                row.names = F, col.names = F, quote = F, sep = "\t")
    write.table(warn, file = paste0(essPrefix, "_WARNINGS_ess.txt"),
                row.names = names(warn), col.names = F, quote = F, sep = "\t", append = TRUE)

  },

  makeMCMCplots = function(mcmc_trace_df, outfile){

    ## all params excluding library alpha params

    pdf(file = outfile)

    opar <- par("mai")

    par("mai" = 0.5 * opar)

    if(shared_active_variance & length(grep(pattern = "^active_variance", colnames(mcmc_trace_df))) > 1)
      mcmc_trace_df <- mcmc_trace_df[,-which(grepl(pattern = "^active_variance", colnames(mcmc_trace_df)))[-1]]

    nparams <- ncol(mcmc_trace_df)

    # if(nparams > 50) stop("too many parameters for plotting")
    if(nparams > 25){

      cat("too many parameters for plotting. Plotting first 25 parameters...\n")
      nparams <- 25

    }

    plot_rows <- ceiling(sqrt(nparams))
    plot_cols <- nparams/plot_rows - (nparams/plot_rows)%%1 + ((nparams/plot_rows)%%1 > 0)

    layout(matrix(seq(plot_cols * plot_rows), ncol = plot_cols, nrow = plot_rows, byrow = T))

    # loop through params and plot
    sapply(seq(ncol(mcmc_trace_df)), function(idx) plot(mcmc_trace_df[,idx], type = "l",
                                          main = colnames(mcmc_trace_df)[idx], cex.main = 0.75,
                                          ylab = "", xlab = ""))

    dev.off()

  },

  makePosteriorPriorPlots = function(mcmc_df, outfile){

    pdf(file = outfile)
    opar <- par("mai")
    par("mai" = 0.5 * opar)
    nparams <- ncol(mcmc_df) - num_libraries + 1 -
      (shared_active_variance) * (num_active_components - 1) -
      1 - (num_active_components == 1)
    plot_rows <- ceiling(sqrt(nparams))
    plot_cols <- nparams/plot_rows - (nparams/plot_rows)%%1 + ((nparams/plot_rows)%%1 > 0)
    layout(matrix(seq(plot_cols * plot_rows), ncol = plot_cols, nrow = plot_rows, byrow = T))

    # mixture means
    post = mcmc_df[,grepl(pattern = "inactive_mean", x = colnames(mcmc_df))]
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    post_density = density(post, from = range[1], to = range[2], cut = 0)
    prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
    prior_density = exp(.self$computeInactiveMeansPriorProbability(prior_xx))
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range,
         xlab = "", ylab = "", main = "inactive_mean", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    post = mcmc_df[,grepl(pattern = "^active_mean", x = colnames(mcmc_df)), drop = F]
    for(k in seq(num_active_components)){
      ifelse(length(threshold_a) == num_active_components, tt <- threshold_a[k], tt <- threshold_a)
      range = c(min(post[,k]) - (max(post[,k]) - min(post[,k])),
                max(post[,k]) + (max(post[,k]) - min(post[,k])))
      post_density = density(post[,k], from = range[1], to = range[2], cut = 0)
      prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
      prior_density = dgamma(prior_xx - tt, active_means_dif_prior_shape, active_means_dif_prior_rate)
      # scale prior
      prior_sf = max(post_density$y) * 0.25 / max(prior_density)
      plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
           main = colnames(post)[k], cex.main = 0.75)
      lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)
    }

    # mixture variances
    post = mcmc_df[,grepl(pattern = "inactive_variance", x = colnames(mcmc_df))]
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    post_density = density(post, from = range[1], to = range[2], cut = 0)

    prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
    prior_density = 1/(prior_xx * (inactive_variances_prior_log_max - inactive_variances_prior_log_min))
    prior_density[prior_xx < 10^inactive_variances_prior_log_min | prior_xx > 10^inactive_variances_prior_log_max] <- 0
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
         main = "inactive_variance", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    post = mcmc_df[,grepl(pattern = "^active_variance", x = colnames(mcmc_df)), drop = F]
    free_variables = 1 + (!shared_active_variance) * (num_active_components - 1)
    for(k in seq(free_variables)){
      range = c(min(post[,k]) - (max(post[,k]) - min(post[,k])),
                max(post[,k]) + (max(post[,k]) - min(post[,k])))
      post_density = density(post[,k], from = range[1], to = range[2], cut = 0)
      prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
      prior_density = 1/((active_variances_prior_log_max - active_variances_prior_log_min) * prior_xx)
      prior_density[prior_xx < 10^active_variances_prior_log_min | prior_xx >10^active_variances_prior_log_max] <- 0
      # scale prior
      prior_sf = max(post_density$y) * 0.25 / max(prior_density)
      plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
           main = colnames(post)[k], cex.main = 0.75)
      lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)
    }

    # mixture weights
    post = mcmc_df[,grepl(pattern = "weight_active", x = colnames(mcmc_df))]
    post_density = density(post, from = 0, to = 1, cut = 0)
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    if(range[1] < 0) range[1] = 0
    if(range[2] > 1) range[2] = 1
    prior_xx = seq(range[1], range[2], by = abs(diff(range))/1000)
    prior_density = exp(.self$computeActiveWeightPriorProbability(prior_xx))
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
         main = "weight_active", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    if(num_active_components > 1){
      post = mcmc_df[,grepl(pattern = "weight_within_active", x = colnames(mcmc_df)), drop = F]
      for(k in seq(num_active_components)){

        ifelse(length(threshold_a) == num_active_components, tt <- threshold_a[k], tt <- threshold_a)

        post_density = density(post[,k], from = 0, to = 1, cut = 0)
        range = c(min(post[,k]) - (max(post[,k]) - min(post[,k])),
                  max(post[,k]) + (max(post[,k]) - min(post[,k])))
        if(range[1] < 0) range[1] = 10^-6
        if(range[2] > 1) range[2] = 1
        prior_xx = seq(range[1], range[2], by = abs(diff(range))/1000)
        prior_density = dbeta(prior_xx, weight_within_active_alpha,
                              weight_within_active_alpha * (num_active_components - 1))
        # scale prior
        prior_sf = max(post_density$y) * 0.25 / max(prior_density)
        plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
             main = colnames(post)[k], cex.main = 0.75)
        lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)
      }
    }


    # spike prob
    post = mcmc_df[,grepl(pattern = "spike", x = colnames(mcmc_df))]
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    if(range[1] < 0) range[1] = 0
    if(range[2] > 1) range[2] = 1
    post_density = density(post, from = range[1], to = range[2], cut = 0)
    prior_xx = seq(range[1], range[2], by = abs(diff(range))/1000)
    prior_density = dbeta(prior_xx, shape1=spike_prior_shape_1, shape2=spike_prior_shape_2)
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
         main = "spike_probability", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    # s0
    post = mcmc_df[,grepl(pattern = "s0", x = colnames(mcmc_df))]
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    post_density = density(post, from = range[1], to = range[2], cut = 0)
    prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
    prior_density = exp(.self$computeS0PriorProbability(prior_xx))
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
         main = "s0", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    # s1
    post = mcmc_df[,grepl(pattern = "s1", x = colnames(mcmc_df))]
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    post_density = density(post, from = range[1], to = range[2], cut = 0)
    prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
    prior_density = exp(.self$computeS1PriorProbability(prior_xx))
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
         main = "s1", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    # tau
    post = mcmc_df[,grepl(pattern = "tau", x = colnames(mcmc_df))]
    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))
    post_density = density(post, from = range[1], to = range[2], cut = 0)
    prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
    prior_density = exp(.self$computeTauPriorProbability(prior_xx))
    # scale prior
    prior_sf = max(post_density$y) * 0.25 / max(prior_density)
    plot(post_density$x, post_density$y, type = "l", xlim = range, xlab = "", ylab = "",
         main = "tau", cex.main = 0.75)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    # alpha_r's
    post = mcmc_df[,grepl(pattern = "alpha_r", x = colnames(mcmc_df)), drop = F]
    density_list_x = vector(mode = "list", length = num_libraries)
    density_list_y = vector(mode = "list", length = num_libraries)
    for(k in seq(num_libraries)){
      ddd <- density(post[,k])
      density_list_x[[k]] <- ddd$x
      density_list_y[[k]] <- ddd$y


    }

    range = c(min(post) - (max(post) - min(post)),
              max(post) + (max(post) - min(post)))

    maxDensity_amongLibs <- max(unlist(density_list_y))

    plot(density_list_x[[1]], density_list_y[[1]], type = "l",
         xlim = range, ylim = c(0, maxDensity_amongLibs),
         xlab = "", ylab = "", main = "alpha_r", cex.main = 0.75)

    for(k in c(2:num_libraries))
      lines(density_list_x[[k]], density_list_y[[k]])

    prior_xx = seq(range[1], range[2], by = abs(diff(range))/100)
    prior_density = dgamma(prior_xx, alpha_r_shape, alpha_r_rate)
    # scale prior
    prior_sf = max(maxDensity_amongLibs) * 0.25 / max(prior_density)
    lines(prior_xx, prior_density * prior_sf, type = "l", lty = 2)


    par("mai" = opar)
    dev.off()
  }

)
