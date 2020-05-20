zigzag$methods(

  findthresholds = function(dat, K = "auto"){

    #get peaks and shoulders
    ps <- .self$get_curve_fitting_list(dat, s = 1)

    if(any(is.na(ps$s))){ #need to fix this terrible code when you have the chance
      ps$s_idx <- c()
      ps$s <- c()
    }

    if(!is.numeric(K)) K <- min(2, sum(ps$s > 0, ps$p > 0))

    #sort shoulders and peaks by the amount of local curvature in the kernel density
    d2y_atshoulders <- ps$d2y[ps$s_idx]
    d2y_atpeaks <- ps$d2y[ps$p_idx]
    shoulders_sorted_by_d2y <- ps$s[order(d2y_atshoulders)]
    peaks_sorted_by_d2y <- ps$p[order(d2y_atpeaks)]

    # get noise peak: max(d2y at peaks/shoulders < 0)
    inactive <- .self$getInactiveThresh(peaks_sorted_by_d2y, shoulders_sorted_by_d2y)

    #get all right peaks up to K number of peaks at log expression > max{0, inactive peak}
    active <- .self$getActiveThresh(peaks_sorted_by_d2y, shoulders_sorted_by_d2y, K)

    allcenters <- sort(c(inactive, active))
    threshes <- allcenters[-length(allcenters)] + diff(allcenters)/2

    return(list(K, threshes))

  },

  getInactiveThresh = function(sorted_peaks_byd2y, sorted_shoulders_byd2y){

    # get noise peak: max(d2y at peaks/shoulders < 0)

    inactive <- c()

    if(sum(sorted_peaks_byd2y < 0) > 0){

      inactive <- sorted_peaks_byd2y[sorted_peaks_byd2y < 0][1]

    }else if(sum(sorted_shoulders_byd2y < 0) > 0){

      inactive <- sorted_shoulders_byd2y[sorted_shoulders_byd2y < 0][1]

    }

    if(length(inactive) == 0) inactive <- 0

    return(inactive)

  },

  getActiveThresh = function(sorted_peaks_byd2y, sorted_shoulders_byd2y, K){

    #get all right peaks up to K number of peaks at log expression > 0

    active <- c()
    num_peak_shoulder_gt0 <- sum(c(sorted_peaks_byd2y, sorted_shoulders_byd2y) > 0)
    num_peak_gt0 <- sum(sorted_peaks_byd2y > 0)

    if(num_peak_shoulder_gt0 >= K){

      if(num_peak_gt0 >= K){

        active <- sorted_peaks_byd2y[sorted_peaks_byd2y > 0][1:K]

      }else{

        active <- c(sorted_peaks_byd2y[sorted_peaks_byd2y > 0],
                    sorted_shoulders_byd2y[sorted_shoulders_byd2y > 0][1:(K - num_peak_gt0)])

      }

    }else{ #if num_peak_shoulder_gt0 < K

      sorted_allcenters <- sort(c(sorted_peaks_byd2y[sorted_peaks_byd2y > 0],
                                  sorted_shoulders_byd2y[sorted_shoulders_byd2y > 0]))

      active <- c(sorted_allcenters,
                  rep(sorted_allcenters[sorted_allcenters == max(sorted_allcenters)],
                      K - length(sorted_allcenters)))

    }

    return(active)

  },

  collapse = function(idx, x){

    x <- sort(x)

    ou <- c(Inf, diff(x))

    i <- which(ou > 1)

    idx[i]

  },

  get_curve_fitting_list = function(x, s=1, lower = min(x), upper = max(x)){

    datrange <- c(sort(x)[0.01 * length(x)], sort(x)[length(x) - 5])

    dens <- density(x, adjust = s, from = lower, to = upper, cut = 0)

    dy <- c(0, diff(dens$y))

    smooth_dy <- predict(loess(dy~dens$x, span = 0.1, normalize = T))

    peaks_idx <- which(diff( smooth_dy >= 0 ) < 0)

    peaks_idx <- peaks_idx[c(dens$x[peaks_idx] > datrange[1] &
                               dens$x[peaks_idx] < datrange[2])]

    d2y <- c(0, diff(dy))

    smooth_d2y <- predict(loess(d2y~dens$x, span = 0.2, normalize = T))

    shoulders_idx <- which(c(0, diff(c(0, diff(smooth_d2y)) > 0 )) > 0)

    shoulders_idx <- shoulders_idx[c(dens$x[shoulders_idx] > datrange[1] &
                                       dens$x[shoulders_idx] < datrange[2])]

    # remove shoulders near peaks
    pairwise_shoulder_peak_dist <- outer(dens$x[shoulders_idx], dens$x[peaks_idx],
                                         FUN = function(X, Y){abs(Y - X)})
    shoulders_idx <- shoulders_idx[rowMins(pairwise_shoulder_peak_dist) > 1]

    # collapse nearby shoulders and peaks to the lower value
    shoulders_idx <- .self$collapse(shoulders_idx, dens$x[shoulders_idx])

    peaks_idx <- .self$collapse(peaks_idx, dens$x[peaks_idx])

    return(list(x = dens$x, y = dens$y, dy = smooth_dy, d2y = smooth_d2y,
                s_idx = shoulders_idx, p_idx = peaks_idx,
                s = dens$x[shoulders_idx], p = dens$x[peaks_idx]))

  }

)
