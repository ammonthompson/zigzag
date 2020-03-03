zigzag$methods(

  drawLibraryDensities = function(x, lib, bl = 0, sf = 1, alpha = 1){

    myx = x

    lines(myx, bl + sf * (1 - weight_active) * (1 - spike_probability) * 1/(1 - (1 - weight_active) * spike_probability) *
            dnorm(myx, mean = inactive_means, sd = sqrt(inactive_variances)),col=rgb(0,0,1,alpha), lwd = 2)

    summed_density <- (1 - weight_active)*(1 - spike_probability)* 1/(1 - (1 - weight_active) * spike_probability) *
      dnorm(myx, mean = inactive_means, sd = sqrt(inactive_variances))

    for(k in 1:num_active_components){

      lines(myx, bl + weight_active * weight_within_active[k] * 1/(1 - (1 - weight_active) * spike_probability) *
              dnorm(myx, mean = active_means[k], sd = sqrt(active_variances[k])), col=rgb(1,(1-1/k),0,alpha), lwd = 2)

      summed_density <- summed_density + weight_active*weight_within_active[k] * 1/(1 - (1 - weight_active) * spike_probability) *
        dnorm(myx, mean = calculate_active_means(active_means_dif)[k], sd=sqrt(active_variances[k]))

    }

    lines(myx, bl + summed_density, col=rgb(0, 1, 0, 0.5 * alpha),lwd = 2)

  },


  drawLibraryHistograms = function(histmids, histdensity, bl = 0, sf = 1){

    histdensity = histdensity * sf

    num_polygons = length(histmids)

    binwidth = histmids[2] - histmids[1]

    for(i in 1:num_polygons){

      rect(xleft = histmids[i]-binwidth/2, xright = histmids[i]+binwidth/2, ybottom = bl, ytop = bl + histdensity[i], col=rgb(0,0,0,0.5))

    }

  },


  plotHistoDensities = function(xy, plotDensity = TRUE){

    xy <- as.matrix(xy)

    num_libs <- ncol(xy)

    histogram_list = list()

    if(length(out_spike_idx) > 2 ){

      for(lib in seq(num_libs)){

         histogram_list[[lib]] <- hist(xy[out_spike_idx, lib], breaks = 100, plot = FALSE)

      }

      max_heights = unlist(lapply(histogram_list, function(libhist){

        return(max(libhist$density[-1]))

      }))

      # print(max_heights)

      total_height = sum(max_heights) - 0.01 * (num_libs - 1)
      max_height = max(max_heights)


      if(num_libs > 1){
        baselines = rev((seq(num_libs) - 1) * (total_height - max_height)/(num_libs - 1))
      }else{
        baselines = c(0)
      }

      plot(NULL, xlim=c(-10,10), ylim = c(0, 1.01 * total_height), xlab = "Y", ylab = "density")

      for(lib in seq(num_libs)){

        .self$drawLibraryHistograms(histogram_list[[lib]]$mids, histogram_list[[lib]]$density, bl = baselines[lib])

        if(plotDensity == TRUE) .self$drawLibraryDensities(histogram_list[[lib]]$mids, lib, bl = baselines[lib])

      }

    }else{
      hist(exp(Yg))
    }

  },


  plotDensities = function(){

    myx=seq(-10,10,by=0.1)

    plot(myx, weight_active*weight_within_active[1]*dnorm(myx, mean = calculate_active_means(active_means_dif)[1],
                                                          sd=sqrt(active_variances[[1]][1])), col="red", type="l",ylim=c(0,0.25),ylab="density",xlab="log TPM")

    for(j in 1:num_libraries){

      lines(density(Xg[,j], bw=0.3), xlim=c(min(myx),max(myx)))

      lines(myx, (1 - weight_active)*(1 - spike_probability)*dnorm(myx, mean = inactive_means, sd = sqrt(inactive_variances)),col="blue", lwd=2)

      summed_density <- (1 - weight_active) * (1 - spike_probability) * dnorm(myx, mean = inactive_means, sd = sqrt(inactive_variances))

      for(k in 1:num_active_components){

        lines(myx, weight_active * weight_within_active[k]*dnorm(myx, mean = calculate_active_means(active_means_dif)[k], sd=sqrt(active_variances[k])), col="red", lwd=2)

        summed_density <- summed_density + weight_active * weight_within_active[k]*dnorm(myx, mean = calculate_active_means(active_means_dif)[k], sd=sqrt(active_variances[k]))

      }

      lines(myx, summed_density, col=rgb(0,1,0,0.3),lwd = 2)

    }

  },


  plotScatterplots = function(){

    plot(Yg, Yg[[2]], xlim=c(-12,10), ylim=c(-12,10), pch =16, col=rgb(0,0,0,0.3))
    points(inactive_means,inactive_means[2], pch = 4, lwd =2, col="blue")
    symbols(inactive_means,inactive_means[2], circles = 1.96 * sqrt(inactive_variances), fg = "blue", inches = F, add = T)
    symbols(inactive_means,inactive_means[2], circles = 1.96 * sqrt(inactive_variances[[2]]), fg = "blue", inches = F, add = T)

    am = calculate_active_means(active_means_dif)
    sapply(1:num_active_components, function(cc){
      points(am[[1]][cc], am[[2]][cc], pch = 4, lwd =2, col="red")
      symbols(am[[1]][cc],am[[2]][cc], circles = (1.96 * sqrt(active_variances[[1]][cc])), fg = "red", inches = F, add = T)
      symbols(am[[1]][cc],am[[2]][cc], circles = (1.96 * sqrt(active_variances[[2]][cc])), fg = "red", inches = F, add = T)
    })



  },


  burninProgressPlot = function(){

    #plot variance trend
    plotdx = seq(-6,8,by=0.1)
    if(variance_g_upper_bound < Inf){
      plot(Yg,variance_g, xlim=c(-6,8), ylim = c(0, variance_g_upper_bound),col=rgb(0,0,0,0.1))
    }else{
      plot(Yg,variance_g, xlim=c(-6,8), ylim = c(0, 0.5 * max(variance_g)),col=rgb(0,0,0,0.1))
    }
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

)
