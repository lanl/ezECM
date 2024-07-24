## © 2024. Triad National Security, LLC. All rights reserved.
##
## This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
##
## (End of Notice)
##
#' Simulate p-values from earthquakes and explosions
#'
#'  p-values are simulated for depth and first polarity discriminants to use in testing classification models.
#'
#'  Methods are adapted from the discriminants listed in \insertCite{anderson2007mathematical}{ezECM}.
#'
#'  # Depth Discriminant
#'
#'  The equation below is used to model p-wave arrival time \eqn{t_i} at seismometer \eqn{i}.
#'
#'  \deqn{t_i = t_0 + T(S_i, S_0) + \epsilon_i}
#'
#'  Where \eqn{t_0} is the time of the event, \eqn{T()} is a function modeling the arrival time (in this case [time_fn]), \eqn{S_i} is the location of seismometer \eqn{i}, \eqn{S_0} is the location of the event, and \eqn{\epsilon_i} is normally distributed error with known variance \eqn{\sigma^2_i}.  Given `N` seismometers, the MLE of the event time \eqn{\hat{t}_0} can be solved as the following:
#'
#'  \deqn{\hat{t}_0 = \frac{\mathrm{tr}\left(\Sigma^{-1} T_i\right) - \mathrm{tr}\left(\Sigma^{-1}T(S,S_0)\right)}{\mathrm{tr}(\Sigma^{-1})}}
#'
#'  Where \eqn{\mathrm{tr}()} is the matrix trace operation, \eqn{\Sigma} is a matrix containing the elements of each \eqn{\sigma_i^2} on the diagonal, \eqn{T_i} is a diagonal matrix containing each \eqn{t_i}, and \eqn{T(S, S_0)} is a diagonal matrix containing each \eqn{T(S_i, S_0)}.  This result is then plugged back into the first equation, which is then used in a normal likelihood function.  Derivatives are taken of the likelihood so that a fast gradient based approch can be used to find the maximum likelihood estimate (MLE) of \eqn{S_0}.
#'
#'  The remainder of the calculation of the p-value is consistent with the *Depth from P-Wave Arrival Times* section of \insertCite{anderson2007mathematical}{ezECM}.  First note \eqn{S_0} is equal to the 3-vector \eqn{(X_0, Y_0, Z_0)^{\top}}.  Given the null hypothesis for the depth of \eqn{\mathrm{H}_0: Z_0 \leq z_0}, the MLE \eqn{(\hat{X}_0, \hat{Y}_0)} given \eqn{Z_0 = z_0} is found.  The sum of squared errors (SSE) is calculated as follows:
#'
#'  \deqn{\mathrm{SSE}(S_0, t_0) = \sum_{i = 1}^N\left(\frac{t_i - t_0 - T(S_i, S_0)}{\sigma_i}\right)^2}
#'
#'  If \eqn{\mathrm{H}_0} is true then the following statistic has a central \eqn{F} distribution with \eqn{1} and \eqn{N-4} degrees of freedom.
#'
#'  \deqn{F_{1,N-4} = \frac{\mathrm{SSE}(S_0, t_0|Z_0 = z_0) - \mathrm{SSE}(S_0, t_0)}{\mathrm{SSE}(S_0, t_0)}}
#'
#'  Because the test has directionality, the \eqn{F} statistic is then converted to a \eqn{T} statistic as such:
#'
#'  \deqn{T_{N-4} = \mathrm{sign}(\hat{Z}_0 - z_0)\sqrt{F_{1,n-4}}}
#'
#'  This \eqn{T} statistic is then used to compute the p-value
#'
#'  #  Polarity of First Motion
#'
#'  Under the null hypothesis that the event is an explosion, (and therefore the true polarity of first motion is always one), the error rate for mistakenly identifying the polarity of first motion is stipulated as the argument `first.polarity$read.err`.  For an error rate of \eqn{\theta} the p-value can then be calculated as follows:
#'
#'  \deqn{\sum_{i = 0}^n {N \choose i} \theta^i(1-\theta)^{N-i}}
#'
#'  Where \eqn{n} is the number of stations where a positive first motion was observed, and \eqn{N} is the total number of stations.
#'
#' @param sims numeric stipulating the number of individual events to generate.
#' @param grid.dim numeric 3-vector providing the extent of the coordinate system in `c(X,Y,Z)` to be used in km.
#' @param seismometer list stipulating variables relating to the seismometers. Providing an incomplete list reverts to the default values.  List elements are:
#'   *  `N` a numeric stipulating the number of seismometers
#'   *  `max.depth` is a numeric providing the maximum depth for a seismometer location in km
#'   *  `Sig` supplying a numeric vector of length `N` to this argument will stipulate the variance in the error in observed arrival time of each of the `N` seismometers.
#'   *  `sig.draws`  a numeric 2-vector, if `Sig` is not provided the variance in arrival time at each station is drawn from [MCMCpack::rinvgamma()] using `shape = sig.draws[1]` and `scale = sig.draws[2]`.
#' @param explosion list stipulating variables regarding a detonation event. Providing an incomplete list reverts to the default values.  List elements are:
#'   *  `max.depth` is a numeric providing the maximum depth of an explosion in km
#'   *  `prob`  is the probability of an explosion, controlling the fraction of events which are explosions.  Value provided must be in the interval \eqn{[0,1]}
#' @param pwave.arrival list stipulating variables regarding the depth from p-wave arrival discriminant. Providing an incomplete list reverts to the default values.  List elements are:
#'   *  `H0` a numeric providing the value of depth in km for the null hypothesis
#'   *  `V`  a numeric stipulating the velocity of p-waves in km.  Used for simulating p-wave arrival times as an argument of [time_fn()].
#'   *  `optim.starts` number of [stats::optim()] starts used to maximize likelihood of event location.
#' @param first.polarity list stipulating variables regarding the depth from the polarity of first motion discriminant.  List elements are:
#'   *  `read.err`  numeric in the interval  \eqn{[0,1]} providing the probability of an error in reading the true first polarity.
#' @returns A data frame, with the number of rows equal to `sims`.  Columns contain the p-values observed for each simulation along with the true event type.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#'   \insertAllCited{}
#'
#' @import MCMCpack
#' @import stats
#'
#' @export
#'
#' @examples
#'
#' test.data <- pval_gen(sims = 5)
#'
#'
pval_gen <- function(sims = 100, grid.dim = c(800, 800, 30),
                     seismometer = list(N = 100, max.depth = 2, Sig = NULL, sig.draws = c(15,2)),
                     explosion = list(max.depth = 5, prob = 0.4),
                     pwave.arrival = list(H0 = 5, V = 5.9, optim.starts = 15),
                     first.polarity  = list(read.err = 0.95)){

  ## Unpacking arguments

  if(any(!is.list(seismometer), !is.list(explosion), !is.list(pwave.arrival), !is.list(first.polarity))){
    warning("Ensure that arguments 'seismometer', 'explosion','pwave.arrival', and 'first.polarity' are provided as a list type data structure.",
            immediate. = TRUE)
  }


  ### Seismometer

  if(is.null(seismometer$N)){
    N <- eval(formals(pval_gen)$seismometer$N)
  }else{
    N <- seismometer$N
  }

  max.depth <- list()

  if(is.null(seismometer$depth)){
    max.depth$seismometer <- eval(formals(pval_gen)$seismometer$max.depth)
  }else{
    max.depth$seismometer <- seismometer$max.depth
  }

  if(is.null(seismometer$Sig)){
    if(is.null(seismometer$sig.draws)){
      Sig <- MCMCpack::rinvgamma(n = N,
                                 shape = eval(formals(pval_gen)$seismometer$sig.draws)[1],
                                 scale = eval(formals(pval_gen)$seismometer$sig.draws)[2])
    }else{
      Sig <- MCMCpack::rinvgamma(n = N,
                                 shape = seismometer$sig.draws[1],
                                 scale = seismometer$sig.draws[2])
    }
  }else{
    if(is.null(dim(seismometer$Sig)) & length(seismometer$Sig) == N){
     Sig <- seismometer$Sig
    }else{
      warning("The specified seismometer$Sig must be a vector of length seismometer$N", immediate. = TRUE)
    }
  }

  ### Explosion

  if(is.null(explosion$max.depth)){
    max.depth$explosion <- eval(formals(pval_gen)$explosion$max.depth)
  }else{
    max.depth$explosion <- explosion$max.depth
  }

  if(is.null(explosion$prob)){
    event.probs <- c(0,1) + eval(formals(pval_gen)$explosion$prob) * c(1,-1)
  }else{
    event.probs <- c(0,1) + explosion$prob * c(1,-1)
  }

  ### pvalue from p-wave arrival time

  if(is.null(pwave.arrival$H0)){
    Z0.cond <- eval(formals(pval_gen)$pwave.arrival$H0)
  }else{
    Z0.cond <- pwave.arrival$H0
  }

  if(is.null(pwave.arrival$V)){
    V.save <- eval(formals(time_fn)$V)
    formals(time_fn)$V <- eval(formals(pval_gen)$pwave.arrival$V)
  }else{
    V.save <- eval(formals(time_fn)$V)
    formals(time_fn)$V <- pwave.arrival$V
  }

  if(is.null(pwave.arrival$optim.starts)){
    starts <- eval(formals(pval_gen)$pwave.arrival$optim.starts)
  }else{
    starts <- pwave.arrival$optim.starts
  }

  ### pvalue from first polarity of motion

  if(is.null(first.polarity$error)){
    read.err <- eval(formals(pval_gen)$first.polarity$read.err)
  }else{
    read.err <- first.polarity$read.err
  }

  eps = sqrt(.Machine$double.eps)

  ### Function for randomly generating sets of labeled p-values ###

  ### This simulation uses a constant number of seismometers in
  ### constant locations with constant variances for varying events.

  p.mat <- data.frame(matrix(NA, ncol = 3, nrow = sims))
  names(p.mat) <- c("p.depth", "p.polarity", "event")
  events <- c("explosion", "earthquake")

  ### Drawing sigmas from inverse gamma distribution


  ### Using random distances from 0-400 in X and Y as well as 0 to 10 in Z

  Si <- matrix(NA, ncol = 3, nrow = N)
  Si[,1] <- grid.dim[1]*stats::runif(N)
  Si[,2] <- grid.dim[2]*stats::runif(N)
  Si[,3] <- max.depth$seismometer*stats::runif(N)

  for(i in 1:sims){

    ### Simulating an event
    #### Sampling event type

    event <- sample(events, size = 1, prob = event.probs)

    #### Generating a random location for the event

    S0 <- stats::runif(1)*grid.dim[1]
    S0 <- c(S0, stats::runif(1)*grid.dim[2])

    if(event == "earthquake"){
      S0 <- c(S0, stats::runif(1)*grid.dim[3])
    }else if(event == "explosion"){
      S0 <- c(S0, stats::runif(1)*max.depth$explosion)
    }


    ##################################################
    ### Depth from p-wave arrival times
    ##################################################

    ## Simulating arrival times

    tn <- P_wave_gen(Si = Si, S0 = S0, Sig = Sig)

    if(is.null(formals(S0_t0_inference)$starts)){
      formals(S0_t0_inference)$starts <- 100
      reset.starts <- TRUE
    }

    ## SSE for the full model
    S0.est.full <- S0_t0_inference(Si = Si,tn = tn, Sig = Sig,
                                   optim.bounds = matrix(c(0,grid.dim[1], 0, grid.dim[2], 0, grid.dim[3]), nrow = 2, ncol = 3),
                                   starts = starts)
    SSE.full <- S0.est.full$SSE

    ## SSE for the model conditioned on a value of Z0 for hypothesis testing
    S0.est.cond <- S0_t0_inference(Si = Si,tn = tn, Sig = Sig, Z0 = Z0.cond,
                                   optim.bounds = matrix(c(0,grid.dim[1], 0, grid.dim[2], 0, grid.dim[3]), nrow = 2, ncol = 3),
                                   starts = starts)

    if(reset.starts){
      formals(S0_t0_inference)$starts <- NULL
    }

    SSE.cond <- S0.est.cond$SSE

    ## Calculating p-value
    f.stat <- (SSE.cond - SSE.full)/(SSE.full/(N-4))
    t.stat <- sign(S0.est.full$est[3] - Z0.cond)*sqrt(f.stat)
    p.depth <- 1-pt(t.stat , df = N - 4)

    ###########################################
    ### Polarity of motion
    ###########################################

    if(event == "explosion"){

      first.polarity <- stats::rbinom(n = N, size = 1, prob = read.err)

      neg.side <- which(first.polarity == 0)
      pos.side <- which(first.polarity == 1)


    }else if(event == "earthquake"){

      ## Specifying a "fault line" along X and Y, which includes the true event location,
      ## to separate the polarities
      ###  Separation of the data using a fault is not currently useful for ECM
      ###  but a model that can use this information may be used in the future

      rndm.pt <- stats::runif(2, min = 0.2, max = 0.8)
      rndm.pt <- rndm.pt*grid.dim[1:2]
      fault.Xy <- rbind(S0[1:2], rndm.pt)
      fault.X <- cbind(1, fault.Xy[,1])
      fault.y <- fault.Xy[,2]
      fault.bhat <- solve(t(fault.X) %*% fault.X) %*% t(fault.X) %*% fault.y

      seis.X <- cbind(1, Si[,1])
      seis.y <- Si[,2]

      yhat <- seis.X %*% fault.bhat

      ## Determining what side of the fault shows positive polarity
      fault.side <- sign(seis.y- yhat)
      fault.side <- sample(c(1,-1), size = 1) * fault.side

      ## Incorporating "reading error"
      fault.side <- fault.side + stats::rbinom(n = N, size = 1, prob = 1-read.err)

      fault.side[fault.side == 2] <- -1
      fault.side[fault.side == 0] <- 1


      neg.side <- which(fault.side == -1)
      pos.side <- which(fault.side == 1)
    }

    n <- length(pos.side)
    p.polarity <- sum(stats::dbinom(0:n, size = N, prob = read.err))

    p.mat[i,] <- c(p.depth,p.polarity, event)

  }

  p.mat[,1:2] <- apply(p.mat[,1:2], 2, as.numeric)


  formals(time_fn)$V <- V.save

  return(p.mat)

}
