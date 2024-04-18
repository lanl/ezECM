## © 2024. Triad National Security, LLC. All rights reserved.
##
## This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
##
## (End of Notice)
##
#' P-wave arrival times simulation
#'
#' Simulates the arrival time of p-waves (without noise) using a mean velocity and euclidian distance, without taking the curvature into account.
#'
#' Used for estimating event location, given a series of seismometer observations.
#'
#' @param X numeric three element vector stipulating the location of a seismometer.  Units of km
#' @param X0 numeric three element vector stipulating the location of an event.  Units of km
#' @param V numeric scalar supplying the p-wave velocity in km/s.  A quick google search (Jan 31st 2023) shows typical values range from 5 to 8 km/s.
#'
#' @return numeric scalar providing the time in seconds for a p-wave to arrive at a seismometer.
#' @export
#'
#' @examples
#'   arrival.time <- time_fn(X = c(100,200,10), X0  = c(500, 80, 25), V = 5)
#'
time_fn <- function(X = NULL, X0 = NULL, V = 7){

  D <- sqrt(sum((X - X0)^2))
  return(D/V)

}

#' Generation of noisy p-wave arrival times
#'
#' Similar utility to [time_fn], however multiple seismometer locations can be provided simultaneously and normally distributed noise is added to the arrival time.
#'
#' @param Si Numeric matrix providing seismometer locations.  Must contain 3 columns corresponding to (X,Y) corrdinates and depth.
#' @param S0 Numeric 3 element vector stipulating the location of an event, elements correspond to (X, Y, Z)
#' @param Sig Numeric vector, or diagonal matrix, providing the variance in observed arrival times at each seismometer.
#' @param neg.obs Logical indicating whether to allow negative observations of time (eg. the observed time of p-wave arrival is before the true time for the event).
#' @param eps Numeric.  If `neg.obs = FALSE` sets of observations are redrawn until all \eqn{t_i - t_0 \leq} `eps`.
#'
#' @return Numeric vector of observation times that correspond to the rows of `Si`
#'
#' @import stats
#'
#' @export
#'
#' @examples
#'
#' pwave.obs <- P_wave_gen(Si = c(100,200,3), S0 = c(400, 500, 4), Sig = 0.05)
#'
P_wave_gen <- function(Si = NULL, S0 = NULL, Sig = NULL, neg.obs = TRUE, eps = sqrt(.Machine$double.eps)){

  if(any(ncol(Si) != 3, is.null(ncol(Si)))){
    if(is.null(dim(Si)) & length(Si) == 3){
      Si <- matrix(Si, ncol = 3, nrow = 1)
    }else{
      warning("Si must be a numeric matrix with 3 columns", immediate. = TRUE)
    }
  }

  if(is.matrix(Sig)){
    if(any(sum(Sig - diag(diag(Sig))) != 0, nrow(Sig) != ncol(Sig), nrow(Sig) != nrow(Si))){
      warning("Sigma must be either a diagonal matrix or a vector of length nrow(Si).")
    }else{
      Sig <- diag(Sig)
    }
  }else{
    if(length(Sig) != nrow(Si)){
      warning("Sigma must be either a diagonal matrix or a vector of length nrow(Si).")
    }
  }

  if(neg.obs == FALSE){
    ti_t0 <- -1
    while(any(ti_t0 <= eps)){
      ti_t0 <- apply(X = Si, 1, time_fn, X0 = S0) + stats::rnorm(nrow(Si), mean = rep(0, times = nrow(Si)), sd = sqrt(Sig))
      }
    }else{
    ti_t0 <- apply(X = Si, 1, time_fn, X0 = S0) + stats::rnorm(nrow(Si), mean = rep(0, times = nrow(Si)), sd = sqrt(Sig))
  }
  return(ti_t0)

}

S0_t0_inference <- function(Si = NULL,tn = NULL, Sig = NULL, Z0 = NULL, starts = NULL,
                            optim.bounds = matrix(c(eps,800, eps, 800, eps, 30), nrow = 2, ncol = 3),
                            eps = sqrt(.Machine$double.eps)){

  params <- list()
  params$Sigi <- 1/Sig
  params$trSigi <- sum(params$Sigi)
  params$Si <- Si
  params$V <- 7
  params$tn <- tn

  if(is.null(Z0)){
    obj <- nll.pwave
    gr <- dnll.pwave
    k <- 3
    outs <- data.frame(matrix(NA, ncol = k + 1, nrow = starts))
    names(outs) <- c("value", "X0", "Y0", "Z0")

    cand <- lhs::randomLHS(n = starts, k = k)
    cand[,1:2] <- cand[,1:2]*max(optim.bounds)
    cand[,3] <- cand[,3]*optim.bounds[2,3]

  }else{
    obj <- nll.pwave.cond
    gr <- dnll.pwave.cond
    k <- 2
    outs <- data.frame(matrix(NA, ncol = k + 1, nrow = starts))
    names(outs) <- c("value", "X0", "Y0")

    params$Z0 <- Z0

    cand <- lhs::randomLHS(n = starts, k = k)
    cand[,1:2] <- cand[,1:2]*max(optim.bounds)
  }


  for(i in 1:starts){

    out <- stats::optim(cand[i,], fn = obj, gr = gr,params = params, method = "L-BFGS-B",
                 lower = optim.bounds[1,1:k], upper = optim.bounds[2,1:k])

    outs$value[i] <- out$value
    outs[i,-1] <- out$par

  }

  best <- unlist(unname(outs[which.min(outs$value), -1]))

  SSE <- 2*obj(best, params = params)

  if(is.null(Z0)){
    T.Si <- apply(Si, 1, time_fn, X0 = best)
  }else{
    T.Si <- apply(Si, 1, time_fn, X0 = c(best, Z0))
  }

  t0hat <- (sum((1/Sig)*tn) -  sum((1/Sig)*T.Si))/sum(1/Sig)

  return(list(est = best, SSE = SSE, t0hat = t0hat))

}


nll.pwave <- function(par, params = NULL){


  Si <- params$Si
  Sigi <- params$Sigi
  tn <- params$tn
  V <- params$V
  trSigi <- params$trSigi

  S0 <- par

  T.Si <- apply(Si, 1, time_fn, X0 = S0)

  t0hat <- (sum(Sigi*tn) - sum(Sigi*T.Si))/trSigi

  p <- (tn - t0hat - T.Si)
  nll <- 0.5*(sum(p^2*Sigi))

  return(nll)
}

dnll.pwave <- function(par, params = NULL){

  Si <- params$Si
  Sigi <- params$Sigi
  tn <- params$tn
  V <- params$V
  trSigi <- params$trSigi

  S0 <- par

  T.Si <- apply(Si, 1, time_fn, X0 = S0)
  Di <- 1/(V*T.Si)
  ViDi <- (1/V*Di)
  ViSigiDi <- ViDi * Sigi

  t0hat <- (sum(Sigi*tn) - sum(Sigi*T.Si))/trSigi
  p <- (tn - t0hat - T.Si)

  dnll <- rep(NA, times = length(S0))

  for(i in 1:length(S0)){

    diff.cord <- S0[i] - Si[,i]

    dnll[i] <- sum((sum(ViSigiDi*diff.cord)/trSigi + ViDi*diff.cord)*Sigi*p)
  }

  return(-dnll)

}

nll.pwave.cond <- function(par, params = NULL){


  Si <- params$Si
  Sigi <- params$Sigi
  tn <- params$tn
  V <- params$V
  trSigi <- params$trSigi
  Z0 <- params$Z0

  S0 <- c(par,Z0)

  T.Si <- apply(Si, 1, time_fn, X0 = S0)

  t0hat <- (sum(Sigi*tn) - sum(Sigi*T.Si))/trSigi

  p <- (tn - t0hat - T.Si)
  nll <- 0.5*(sum(p^2*Sigi))

  return(nll)
}

dnll.pwave.cond <- function(par, params = NULL){

  Si <- params$Si
  Sigi <- params$Sigi
  tn <- params$tn
  V <- params$V
  trSigi <- params$trSigi
  Z0 <- params$Z0

  S0 <- c(par,Z0)

  T.Si <- apply(Si, 1, time_fn, X0 = S0)
  Di <- 1/(V*T.Si)
  ViDi <- (1/V*Di)
  ViSigiDi <- ViDi * Sigi

  t0hat <- (sum(Sigi*tn) - sum(Sigi*T.Si))/trSigi
  p <- (tn - t0hat - T.Si)

  dnll <- rep(NA, times = length(par))

  for(i in 1:length(par)){

    diff.cord <- S0[i] - Si[,i]

    dnll[i] <- sum((sum(ViSigiDi*diff.cord)/trSigi + ViDi*diff.cord)*Sigi*p)
  }

  return(-dnll)

}
