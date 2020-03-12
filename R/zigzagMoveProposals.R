zigzag$methods(

  scale_move = function(param, tuningParam){

    scale_factor <- exp(tuningParam * (runif(length(param)) - 0.5))

    proposed_param <- param * scale_factor

    # return scale_factor for computing hastings ratio

    return(list(proposed_param, scale_factor))

  },

  positive_slide_move = function(param, upper = Inf, tuningParam){

    new_param <- abs(param + rnorm(1, 0, tuningParam))

    if(new_param > upper){

      new_param <- 2 - new_param

    }

    return(new_param)

  },

  two_boundary_slide_move = function(param, lower, upper, tuningParam){

    #new_param <- param + rnorm(1, 0, tuningParam)
    #new_param <- param + runif(1, 0, tuningParam) - tuningParam/2
    new_param <- runif(1, param - tuningParam, param + tuningParam)

    if(new_param > upper){

      new_param <- 2 * upper - new_param

    }

    new_param <- lower + abs(new_param - lower)


    return(new_param)

  },

  no_boundery_slide_move = function(param, tuningParam){

    return(param + rnorm(1, 0, tuningParam))

  },

  twoSigma_slide_move = function(params, lower, upper, tuningParam){

    delta <- rnorm(1, 0, tuningParam)

    new_param1 <- params[1] + delta

    new_param2 <- params[2] - delta

    if(new_param1 > upper){

      new_param1 <- 2 * upper - new_param1

    }

    if(new_param2 > upper){

      new_param2 <- 2 * upper - new_param2

    }

    new_param1 <- lower + abs(new_param1 - lower)

    new_param2 <- lower + abs(new_param2 - lower)

    return(c(new_param1, new_param2))

  },

  twoMean_slide_move = function(params, tuningParam){

    delta <- rnorm(1, 0, tuningParam)

    newparam1 <- abs(params[1] + delta)

    newparam2 <- abs(params[2] - delta)

    return(c(newparam1, newparam2))

  }


)
