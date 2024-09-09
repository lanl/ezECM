## © 2024. Triad National Security, LLC. All rights reserved.
##
## This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
##
## (End of Notice)


#' New Event Categorization With Bayesian Inference
#'
#' @param object an object of `class` `"BayesECM"` obtained as the trained model output using the [BayesECM()] function.
#' @param Ytilde `data.frame` of unlabeled observations to be categorized.  Must contain the same discriminant names as the training data used in the provided "BayesECM" object.  Each row is an individual observation.  Missing data is specified with `NA`
#' @param thinning integer, scalar.  Values greater than one can be provided to reduce computation time.  See details.
#' @param mixture_weights character string describing the weights of the distributions in the mixture to be used for prediction.  The default,`"training"` will utilize weights according to likelihood and prior specifications, while supplying the string `"equal"` will assume the marginal predictive distribution of each category is independent of the data and utilize equal weights.
#' @param ... not used
#'
#' @return Returns a `list`.  The list element `epz` is a matrix with `nrow(Ytilde)` rows, corresponding to each event used for prediction, and `K` named columns.  Each column of `epz` is the expected category probability of the row stipulated event.  The remainder of the list elements hold data including `Ytilde`, information about additonal variables passed to `predict.BayesECM`, and data related to the previous [BayesECM()] fit.
#' @export
#'
#' @details
#'
#' The data in `Ytilde` should be the p-values \eqn{\in (0,1]}.  The transformation applied to the data used to generate `object` is automatically applied to `Ytilde` within the `predict.BayesECM()` function.
#'
#' For a given event with an unknown category, a Bayesian ECM model seeks to predict the expected value of the latent variable \eqn{\tilde{\mathbf{z}}_K}, where \eqn{\tilde{\mathbf{z}}_K} is a vector of the length \eqn{K}, and \eqn{K} is the number of event categories.  A single observation of \eqn{\tilde{\mathbf{z}}_K} is a draw from a \href{https://en.wikipedia.org/wiki/Categorical_distribution}{Categorical Distribution}.
#'
#' The expected probabilities stipulated within the categorical distribution of \eqn{\tilde{\mathbf{z}}_K} are conditioned on any imputed missing data, prior hyperparameters, and individually each row of `Ytilde`. The output from [predict.BayesECM()] are draws from the distribution of \eqn{\mathbf{E}[\tilde{\mathbf{z}}_K|\tilde{\mathbf{y}}_{\tilde{p}}, \mathbf{Y}^{+}, \mathbf{\eta}, \mathbf{\Psi}, \mathbf{\nu}, \mathbf{\alpha}] = p(\tilde{\mathbf{z}}_K|\tilde{\mathbf{y}}_{\tilde{p}}, \mathbf{Y}^{+}, \mathbf{\eta}, \mathbf{\Psi}, \mathbf{\nu}, \mathbf{\alpha})}, where \eqn{\mathbf{Y}^{+}} represents the observed values within the training data.
#'
#' The argument `mixture_weights` controls the value of \eqn{p(\tilde{\mathbf{z}}_K|\mathbf{Y}_{N \times p}, \mathbf{\alpha})}, the probability of each \eqn{\tilde{z}_k = 1}, before \eqn{\tilde{\mathbf{y}}_{\tilde{p}}} is observed.  The standard result is obtained from the prior hyperparameter values in \eqn{\mathbf{\alpha}} and the number of unique events in each \eqn{\mathbf{Y}_{N_k \times p}}.  Setting `mixture_weights = "training"` will utilize this standard result in prediction.  If the frequency of the number events used for each category in training is thought to be problematic, providing the argument `mixture_weights = "equal"` sets \eqn{p(\tilde{z}_1 = 1|\mathbf{Y}_{N \times p}) =  \dots = p(\tilde{z}_K = 1|\mathbf{Y}_{N \times p}) = 1/K}.  If the user wants to use a set of \eqn{p(\tilde{z}_k = 1|\mathbf{Y}_{N \times p})} which are not equal but also not informed by the data, we suggest setting the elements of the hyperparameter vector \eqn{\mathbf{\alpha}} equal to values with a large magnitude and in the desired ratios for each category.  However, this can cause undesirable results in prediction if the magnitude of some elements of \eqn{\mathbf{\alpha}} are orders larger than others.
#'
#' To save computation time, the user can specify an integer value for `thinning` greater than one.  Every `thinning`th Markov-chain Monte-Carlo sample is used for prediction.  This lets the user take a large number of samples during the training step, allowing for better mixing.  See details in a package vignette by running \code{vignette("syn-data-code", package = "ezECM")}
#'
#' @examples
#'
#' csv_use <- "good_training.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' training_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' trained_model <- BayesECM(Y = training_data, BT = c(10,1000))
#'
#' csv_use <- "good_newdata.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' new_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' bayes_pred <- predict(trained_model,  Ytilde = new_data)
#'
#'
#' @method predict BayesECM
#'
predict.BayesECM <- function(object, Ytilde, thinning = 1, mixture_weights = "training", ...){


  Ytilde$event <- NULL
  if(!all(names(Ytilde) %in% names(object$Y[[1]]))){
    missing.names <- names(Ytilde)[!(names(Ytilde) %in% names(object$Y))]
    stop(paste0("Error: the discriminant names ", paste(missing.names, collapse = ", "), " are not found in the training data names of:  ", paste(names(object$Y[[1]]))))
  }

  Ytilde <- Ytilde[names(object$Y[[1]])]


  if(object$data$transform == "logit"){
    Ytilde <- log(Ytilde) - log(1-Ytilde)
  }else if(object$data$transform == "arcsin"){
    Ytilde <- (2/pi) * asin(sqrt(Ytilde))
  }else{
    warning("Transform data missing from the supplied BayesECM object.  Ytilde will not be transformed for prediction.", immediate. = TRUE)
  }


  if(mixture_weights == "equal"){
    eq_wts <- TRUE
  }else if(mixture_weights == "training"){
    eq_wts <- FALSE
  }else{
    stop(paste0("Error: the argument 'mixture_weights' must be either 'training' or 'equal'. ", mixture_weights, " is not a valid specification."))
  }

 epz <- BayesECM_pred(BayesECM_obj = object, Ytilde = Ytilde, thinning = thinning, eq_wts = eq_wts)

 predout <- structure(list(epz = epz$pz, Ytilde = Ytilde, Ytilde_group = epz$Ytilde_group , umissing = epz$umissing, thinning = epz$thinning, BayesECMfit = object), class = "BayesECMpred")
 return(predout)


}


#' Bayesian Event Matrix Categorization
#'
#' Training a Bayesian ECM (B-ECM) model
#'
#' @details
#'
#' The output of `BayesECM()` provides a trained Bayesian Event Categorization Matrix (B-ECM) model, utilizing the data and prior parameter settings .  If there are missing values in `Y`, these values are imputed.  A trained `BayesECM` model is then used with the [predict.BayesECM()] function to calculate expected category probabilities.
#'
#' ##  Data Preparation
#'
#' Before the data in `Y` is used with the model, the p-values \eqn{\in (0,1]} are transformed in an effort to better align the data with some properties of the normal distribution.  When `transform == "logit"` the inverse of the logistic function \eqn{Y_{N \times p} = \log\left(\texttt{Y}\right) - \log\left(1-\texttt{Y}\right)} maps the values to the real number line.  Values of `Y` exactly equal to 0 or 1 cannot be used when `transform == "logit"`.  Setting the argument `transform == "arcsin"` uses the transformation \eqn{Y_{N\times p} = 2/\pi \times \mathrm{arcsin}\sqrt{Y}} further described in \insertCite{anderson2007mathematical;textual}{ezECM}.  From here forward, the variable \eqn{Y_{N \times p}} should be understood to be the transformation of `Y`, where \eqn{N} is the total number of rows in `Y` and \eqn{p} is the number of discriminant columns in `Y`.
#'
#'  ## The Model
#'
#' The B-ECM model structure can be found in a future publication, with some details from this publication are reproduced here.
#'
#' B-ECM assumes that all data is generated using a mixture of \eqn{K} normal distributions, where \eqn{K} is equal to the number of unique event categories. Each component of the mixture has a unique mean of \eqn{\mu_k}, and covariance of \eqn{\Sigma_k}, where \eqn{k \in \{1 , \dots , K\}} indexes the mixture component.  The likelihood of the \eqn{i^{\mathrm{th}}} event observation \eqn{y^i_p} of \eqn{p} discriminants can be written as the sum below.
#'
#' \deqn{\sum_{k = 1}^K \pi_k \mathcal{N}(y^i_p; \mu_k, \Sigma_k)}
#'
#' Each Gaussian distribution in the sum is weighted by the scalar variable \eqn{\pi_k}, where \eqn{\sum_{k=1}^K \pi_k =1} so that the density integrates to 1.
#'
#' There are prior distributions on each \eqn{\mu_k, \Sigma_k}, and \eqn{\pi}, where \eqn{\pi} is the vector of mixture weights \eqn{\{\pi_1, \dots , \pi_K\}}.  These prior distributions are detailed below.  These parameters are important for understanding the model, however they are integrated out analytically to reduce computation time, resulting in a marginal likelihood \eqn{p(Y_{N_k \times p}|\eta_k, \Psi_k, \nu_k)} which is a mixture of matrix t-distributions.  \eqn{Y_{N_k \times p}} is a matrix of the total data for the \eqn{k^{\mathrm{th}}} event category containing \eqn{N_k} total event observations for training.  The totality of the training data can be written as \eqn{Y_{N \times p}}, where \eqn{N = N_1 + \dots + N_K}.
#'
#' `BayesECM()` can handle observations where some of the \eqn{p} discriminants of an observation are missing.  The properties of the conditional matrix t-distribution are used to impute the missing values, thereby accounting for the uncertainty related to the missing data.
#'
#' ##  Prior Distributions
#'
#'  The posterior distributions  \eqn{p(\mu_k|Y_{N_k \times p}, \eta_k)}, \eqn{p(\Sigma_k|Y_{N_k \times p}, \Psi_k, \nu_k)}, and \eqn{p(\pi|Y_{N \times p}, \alpha)} are dependent on the specifications of prior distributions \eqn{p(\mu_k|\Sigma_k, \eta_k)}, \eqn{p(\Sigma_k| \Psi_k, \nu_k)}, and \eqn{p(\pi|\alpha)}.
#'
#'  \eqn{p(\mu_k|\Sigma_k, \eta_k)} is a multivariate normal distribution with a mean vector of \eqn{\eta_k} and is conditional on the covariance \eqn{\Sigma_k}.   \eqn{p(\Sigma_k|\Psi_k, \nu_k)} is an \href{https://en.wikipedia.org/wiki/Wishart_distribution}{Inverse Wishart} distribution with degrees of freedom parameter \eqn{\nu_k}, or `nu`, and scale matrix \eqn{\Psi_k}, or `Psi`.  \eqn{p(\pi|\alpha)} is a \href{https://en.wikipedia.org/wiki/Dirichlet_distribution}{Dirichlet distribution} with the parameter vector \eqn{\alpha} of length \eqn{K}.
#'
#' The ability to use `"default"` priors has been included for ease of use with various settings of the `priors` function argument.  The default prior hyperparameter values differ for the argument of `transform` used, and the values can be inspected by examining the output of the `BayesECM()` function.  Simply setting `priors = "default"` provides the same default values for all \eqn{\eta_k, \Psi_k, \nu_k} in the mixture.  If all prior parameters are to be shared between all event categories, but some non-default values are desirable then supplying a list of a similar structure as `priors = list(eta = rep(0, times = ncol(Y) - 1), Psi = "default", nu = "default", alpha = 10)` can be used, where setting a list element `"default"` can be exchanged for the correct data structure for the relevant data structure.
#'
#'  If one wishes to use some default values, but not share all parameter values between each event category, or wishes to specify each parameter value individually with no defaults, we suggest running and saving the output `BayesECM(Y = Y, BT = c(1,2))$priors`.  Note that when specifying `eta` or `Psi` it is necessary that the row and column order of the supplied values corresponds to the column order of `Y`.
#'
#'
#'
#' @param Y `data.frame` of training data, with rows corresponding to \eqn{N} individual observations, and columns corresponding to \eqn{p} discriminants.  An additional column named `"event"` is required, which labels each row with the known event category.  Missing data is specified with `NA`.  `dim(Y)` must equal `c(N, p + 1)`.
#' @param BT integer vector of length 2, stipulating the number of `c("Burn-in", "Total")` Markov-Chain Monte-Carlo samples drawn.
#' @param priors list of parameters to be used in prior distributions.  See details.
#' @param verb logical.  A setting of `TRUE` prints a progress bar as samples are drawn, as well as warnings when they occur.
#' @param transform character string specifying the transform to use on the elements of `Y` before fitting the model.  Options are `"logit"` and `"arcsin"` with `"logit"` being the default.  See details.
#'
#' @return
#'
#' Returns an object of `class("BayesECM")`.  If there are missing data in the supplied argument `Y` the object contains Markov-chain Monte-Carlo samples of the imputed missing data.  Prior distribution parameters used are always included in the output.  The primary use of an object returned from [BayesECM()] is to later use this object to categorize unlabeled data with the [predict.BayesECM()] function.
#'
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @examples
#'
#' csv_use <- "good_training.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' training_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' trained_model <- BayesECM(Y = training_data)
#'
#'
BayesECM <- function(Y, BT = c(100, 1000), priors = "default", verb = FALSE, transform = "logit"){

  if(!("event" %in% names(Y))){
    stop("Supplied Y data.frame must contain a column named 'event' which specifies the data labels")
  }
  if(!(transform %in% c("logit", "arcsin")) | length(transform) != 1 | !is.character(transform) | is.null(transform)){
    stop("Supplied argument 'transform' must be a single character string equal to 'logit' or 'arcsin'.")
  }
  Ylist <- list()
  cats <- unique(Y$event)

  if(length(cats) <= 1){
    stop("Number of event categories must be 2 or greater.")
  }

  discriminants <- names(Y)[-which(names(Y) == "event")]

  for(i in cats){
    Ylist[[i]] <- unname(as.matrix(Y[Y$event == i, -which(names(Y) == "event")]))
    if(transform == "logit"){
      Ylist[[i]] <- log(Ylist[[i]]) - log(1-Ylist[[i]])
    }else{
      Ylist[[i]] <- (2/pi) * asin(sqrt(Ylist[[i]]))
    }
  }

  cats <- sort(cats)
  Ylist <- Ylist[cats]

  discrim_warn <- sapply(Ylist, function(X){
    apply(X,2,function(X){
      all(is.na(X))
    })
  })

  if(verb){
    if(sum(discrim_warn) != 0){
      for(i in 1:length(cats)){
        if(sum(discrim_warn[,cats[i]]) != 0){
          discrim_missing <- discriminants[which(discrim_warn[,cats[i]])]
          warning(paste0("Discriminants ", paste(discrim_missing, collapse = ", "), " fully missing from event category ", cats[i], "\n", "Inference for this missing data only informed by prior distributions.\n"), immediate. = TRUE)
        }
      }
    }
  }



  if(all(cats %in% names(priors))){
    priors <- priors[cats]
  }

  prior_check_message <- prior_checks(priors = priors, Y = Ylist)

    if(!is.null(prior_check_message)){
      stop(prior_check_message)
    }

  priors <- prior_allocation(priors = priors, Y = Ylist, transform = transform)


  Bayes_train <- BayesECM_train(Y = Ylist, BT = BT, priors = priors, verb = verb, cats = cats)

  for(i in 1:length(Bayes_train$Y)){
    Bayes_train$Y[[i]] <- data.frame(Bayes_train$Y[[i]])
    names(Bayes_train$Y[[i]]) <-  discriminants
  }

  Bayes_train$data$transform <- transform
  return(train = Bayes_train)

}



prior_checks <- function(priors = NULL, Y = NULL){
  if(all(priors == "default")){
    return()
  }

  N <- unname(sapply(Y, nrow))
  p <- max(sapply(Y, ncol))
  K <- length(Y)

  if(!(length(priors) %in% c(4,K))){
    return("The 'priors' argument must be set to the character string 'default', must contain shared values for each event category, or be of length equal to the number of event categories.  See documentation.")
  }

  if(length(priors) == K){
    if(any(unlist(priors) == "default")){
      return("Default priors for any parameters cannot be used if length(priors) is equal to the number of event categories.")
    }
    for(k in 1:K){
      if(!all(names(priors[[k]]) %in% c("eta", "Psi", "alpha", "nu"))){
        return(paste0("Some prior parameters missing.  Check priors[[", k, "]]"))
      }
      if(!isSymmetric(priors[[k]]$Psi)){
        return(paste0("All supplied Psi matrices must be square and symmetric.  Check priors[[", k, "]]$Psi"))
      }else if((nrow(priors[[k]]$Psi) != p)){
        return(paste0("Supplied priors[[", k, "]]$Psi must be of dimension ", p, " for the supplied data set."))
      }else if(!all(eigen(priors[[k]]$Psi)$values > 0)){
        return(paste0("Supplied priors[[", k, "]]$Psi must be full rank and positive definite."))
      }

      priors[[k]]$eta <- as.vector(priors[[k]]$eta)
      if(length(priors[[k]]$eta) != p){
        return(paste0("priors[[", k, "]]$eta must be a vector of length ", p, " for this data set."))
      }

      if(length(priors[[k]]$alpha) != 1 | !is.numeric(priors[[k]]$alpha) | priors[[k]]$alpha <= 0){
        return(paste0("priors[[", k, "]]$alpha must be a single numeric value grater than 0."))
      }

      if(length(priors[[k]]$nu) != 1 | !is.numeric(priors[[k]]$nu) | priors[[k]]$nu < (p-1)){
        return(paste0("priors[[", k, "]]$nu must be a single numeric value grater than ", p-1, " preferably greater than ", p))
      }

    }
  }

  if(length(priors) == 4){
    if(!all(names(priors) %in% c("eta", "Psi", "alpha", "nu"))){
      return("Some prior parameters missing")
    }
    if(!isSymmetric(priors$Psi)){
      return("All supplied Psi matrices must be square and symmetric.")
    }else if((nrow(priors$Psi) != p)){
      return(paste0("Supplied priors$Psi must be of dimension ", p, " for the supplied data set."))
    }else if(any(eigen(priors$Psi)$values <= 0)){
      return("Supplied priors$Psi must be full rank and positive definite.")
    }

    priors$eta <- as.vector(priors$eta)
    if(length(priors$eta) != p){
      return(paste0("priors$eta must be a vector of length ", p, " for this data set."))
    }

    if(length(priors$alpha) != 1 | !is.numeric(priors$alpha) | priors$alpha <= 0){
      return(paste0("priors$alpha must be a single numeric value grater than 0."))
    }

    if(length(priors$nu) != 1 | !is.numeric(priors$nu) | priors$nu < (p-1)){
      return(paste0("priors$nu must be a single numeric value grater than ", p-1, " preferably greater than ", p))
    }

  }
  return()
}


prior_allocation <- function(priors = NULL, Y = NULL, transform = NULL){


  N <- unname(sapply(Y, nrow))
  p <- max(sapply(Y, ncol))
  K <- length(Y)
  eps <- sqrt(.Machine$double.eps)
  discrim_names <- names(Y[[1]])

  if(any(priors == "default")){
    ## Sets default values to parameters specified as default
    if(priors == "default"){
      priors <- list(eta = "default", Psi = "default", nu = "default", alpha = "default")
    }
    default_prior_names <- names(which(priors == "default"))
    if("eta" %in% default_prior_names){
      if(transform == "logit"){
        eta_temp <- rep(0, times = p)
      }else if(transform == "arcsin"){
        eta_temp <- rep(0.5, times = p)
      }else{
        stop("Error on prior allocation with respect to the prior specific to the transform.")
      }
    }else{
      eta_temp <- priors$eta
    }
    names(eta_temp) <- discrim_names

    if("Psi" %in% default_prior_names){
      if(transform == "logit"){
        Psi_temp <- diag(p)
      }else if(transform == "arcsin"){
        Psi_temp <- diag(p)*0.1
      }else{
        stop("Error on prior allocation with respect to the prior specific to the transform.")
      }

    }else{
      Psi_temp <- priors$Psi
    }
    colnames(Psi_temp) <- discrim_names
    rownames(Psi_temp) <- discrim_names

    if("nu" %in% default_prior_names){
      nu_temp <- p
    }else{
      nu_temp <- priors$nu
    }

    if("alpha" %in% default_prior_names){
      alpha_temp <- 1/2
    }else{
      alpha_temp <- priors$alpha
    }

    priors <- list()

    for(k in 1:K){
      priors[[k]] <- list(eta = eta_temp, Psi = Psi_temp, nu = nu_temp, alpha = alpha_temp)
    }

  }else if(all(c("eta", "Psi", "nu", "alpha") %in% names(priors))){

    ## If user specified values, shared between all components of the mixture, are supplied, these values are used.

    names(priors[["eta"]]) <- discrim_names
    rownames(priors[["Psi"]]) <- discrim_names
    colnames(priors[["Psi"]]) <- discrim_names
    priors_temp <- priors[c("eta", "Psi", "nu", "alpha")]
    priors <- list()

    for(k in 1:K){
      priors[[k]] <- priors_temp
    }

  }

  ## If the previous two if statements are not triggered, and no warnings were given from the prior_checks function,
  ## it is assumed that the user supplied the correct format for the prior list, and this is used

  names(priors) <- names(Y)

  return(priors)

}

BayesECM_train <- function(Y = NULL, BT = NULL, priors = NULL, verb = NULL, cats = NULL){


  ## Returns imputed data, some data structures that help with prediction, priors used

  # General Data
  N <- unname(sapply(Y, nrow))
  p <- max(sapply(Y, ncol))
  K <- length(Y)
  pvec <- 1:p
  Ip <- diag(p)
  Yimpute <- list()

  Ydrawn <- Y
  NSighat <- (N + 2)/(N + 1)

  # Save for each event type which columns have missing data, the number of missing data, and the index of missing data

  N_minus <- N_plus <- list()
  Y_minus_index <- list()
  IpJ <- list()
  I_one_N_J <- list()
  # I_N - (1/(N + 1))*J_N
  Im1N1J <- list()
  Psi22iPsi21 <- list()
  Psi11cond <- list()
  Ybar <- list()
  predmean <- SighatN <- list()
  p_part_draws <- p_full_draws <- list()

  for(k in 1:K){
    Psi22iPsi21[[k]] <- list()
    Psi11cond[[k]] <- rep(NA, times = p)
    Y_minus_index[[k]] <- apply(Y[[k]], 2, function(X){which(is.na(X))})
    N_minus[[k]] <- apply(Y[[k]],2, function(X){sum(is.na(X))})
    N_plus[[k]] <- N[k] - N_minus[[k]]
    IpJ[[k]] <- diag(N[k]) + matrix(1, ncol = N[k], nrow = N[k])
    Yimpute[[k]] <- list()
    Ybar[[k]] <- matrix(NA, ncol = BT[2], nrow = p)
    SighatN[[k]] <-  matrix(NA, ncol = BT[2], nrow = p*(1+p)/2)
    predmean[[k]] <- matrix(NA, ncol = BT[2], nrow = p)
    I_one_N_J[[k]] <- diag(N[k]) - (1/(N[k]+ 1))*matrix(1, ncol = N[k], nrow = N[k])
    p_part_draws[[k]] <- which(!(N_plus[[k]] %in% c(0,N[k])))

    if(0 %in% N_plus[[k]]){
      p_full_draws[[k]] <- (1:p)[N_plus[[k]] == 0]
    }else{
      p_full_draws[[k]] <- numeric(0)
    }
    for(j in 1:p){
      A <- Ip[,c(j, pvec[-j])]
      Psitemp <- priors[[k]]$Psi
      Psitemp <- t(A) %*% Psitemp %*% A
      Psi22 <- Psitemp[2:p, 2:p]
      Psi21 <- Psitemp[2:p,1]
      Psi22iPsi21[[k]][[j]] <- Psi22iPsi21temp <- solve(Psi22, Psi21)
      Psi11cond[[k]][j] <- drop(Psitemp[1,1] - t(Psi21) %*% Psi22iPsi21temp)
      Yimpute[[k]][[j]] <- matrix(NA, ncol = BT[2], nrow = N_minus[[k]][j])
      Ydrawn[[k]][is.na(Ydrawn[[k]][,j]),j] <- priors[[k]]$eta[j]
    }
  }


  Kfull <- rep(NA,times = K)
  for(k in 1:K){
    Kfull[k] <- all(N[k] == N_plus[[k]])
  }

  Kmissing <- which(!Kfull)
  Kfull <- which(Kfull)



  if(verb & length(Kmissing) != 0){
    print(paste0("Training, gathering ", BT[2], " Monte-Carlo samples."))
    pb <- utils::txtProgressBar(min = 1, max = BT[2], initial = 1, style = 3)
    pbupdate <- floor(seq(from = 0, to = BT[2], length.out = 100))
    on.exit(close(pb))
  }

  for(i in 1:BT[2]){

    for(k in Kmissing){

      # Drawing discriminants for which there is partial data
      for(j in p_part_draws[[k]]){



        y2 <- Ydrawn[[k]][-Y_minus_index[[k]][[j]], j, drop = FALSE]
        y3 <- Ydrawn[[k]][Y_minus_index[[k]][[j]], -j, drop = FALSE]
        y4 <- Ydrawn[[k]][-Y_minus_index[[k]][[j]], -j, drop = FALSE]


        y1 <- drop(rconditional_Mt(y2 = y2, y3 = y3, y4 = y4, nu = priors[[k]]$nu, Np = N_plus[[k]][j],
                              etap = priors[[k]]$eta[-j], etam = priors[[k]]$eta[j], Psi22 = priors[[k]]$Psi[-j, -j],
                              Psi22iPsi21 = Psi22iPsi21[[k]][[j]], Omega = Psi11cond[[k]][j]))


        Ydrawn[[k]][Y_minus_index[[k]][[j]], j] <- y1
        Yimpute[[k]][[j]][,i] <- y1
      }
      # Drawing discriminants for which all data is missing
      for(j in p_full_draws[[k]]){

        y3y4 <- Ydrawn[[k]][, -j, drop = FALSE]

        y1y2 <- drop(rcol_conditional_Mt(y3y4 = y3y4, nu = priors[[k]]$nu, etap = priors[[k]]$eta[-j],
                                    etam = priors[[k]]$eta[j], Psi22 = priors[[k]]$Psi[-j, -j],
                                    Psi22iPsi21 = Psi22iPsi21[[k]][[j]], Omega = Psi11cond[[k]][j]))


        Ydrawn[[k]][,j] <- y1y2
        Yimpute[[k]][[j]][,i] <- y1y2
      }

      Ybar[[k]][,i] <- apply(Ydrawn[[k]],2, mean)

      Yswept <- sweep(Ydrawn[[k]], 2, priors[[k]]$eta)


      SighatNdraw_temp <- (priors[[k]]$Psi + t(Yswept) %*% I_one_N_J[[k]] %*% Yswept) * NSighat[k]
      SighatNdraw_temp <- (SighatNdraw_temp + t(SighatNdraw_temp))/2
      SighatN[[k]][,i] <- SighatNdraw_temp[upper.tri(SighatNdraw_temp, diag = TRUE)]

      predmean[[k]][,i] <- (1/(N[k] + 1))*(N[k]*Ybar[[k]][,i] + priors[[k]]$eta)

    }

    if(verb & length(Kmissing) != 0){
      if(i %in% pbupdate){
        utils::setTxtProgressBar(pb, i)
      }
    }
  }

  ## Remove Burn in


    for(k in Kmissing){
      for(j in 1:p){
      Yimpute[[k]][[j]] <- Yimpute[[k]][[j]][,-(1:BT[1])]
      }
      Ybar[[k]] <- Ybar[[k]][,-(1:BT[1])]
      SighatN[[k]] <- SighatN[[k]][, -(1:BT[1])]
      predmean[[k]] <- predmean[[k]][, -(1:BT[1])]
    }

  for(k in Kfull){
    Ybar[[k]] <- apply(Ydrawn[[k]],2, mean)

    Yswept <- sweep(Ydrawn[[k]], 2, priors[[k]]$eta)

    SighatNdraw_temp <- (priors[[k]]$Psi + t(Yswept) %*% I_one_N_J[[k]] %*% Yswept) * NSighat[k]
    SighatNdraw_temp <- (SighatNdraw_temp + t(SighatNdraw_temp))/2
    SighatN[[k]] <- SighatNdraw_temp[upper.tri(SighatNdraw_temp, diag = TRUE)]

    predmean[[k]] <- (1/(N[k] + 1))*(N[k]*Ybar[[k]] + priors[[k]]$eta)
  }

  names(predmean) <- names(Y)

  alpha <- sapply(priors, function(X){X$alpha})

  a_p <- alpha + N

  return(structure(list(Y = Y,
                        MCMC = list(Ybar = Ybar, SighatN = SighatN, Yimpute = Yimpute, la_p = log(a_p), a_p  = a_p, predmean = predmean),
                        priors = priors,
                        data = list(Y_minus_index = Y_minus_index, Kfull = Kfull)),
                   class = "BayesECM"))


}

rconditional_Mt <- function(y2 = NULL, y3 = NULL, y4 = NULL, nu = NULL, Np = NULL,
                            etap = NULL, etam = NULL, Psi22 = NULL, Psi22iPsi21 = NULL, Omega = NULL){

  y3_swept <- sweep(y3, 2, etap)
  y4_swept <- sweep(y4, 2, etap)

  Psi22ity3_swept <- solve(Psi22, t(y3_swept))
  Psi22ity4_swept <- solve(Psi22, t(y4_swept))

  M1r <- etam + y3_swept %*% Psi22iPsi21
  M2r <- etam + y4_swept %*% Psi22iPsi21
  Sig12 <- y3_swept %*% Psi22ity4_swept + 1
  Sig22 <- y4_swept %*% Psi22ity4_swept + 1
  diag(Sig22) <- diag(Sig22) + 1
  Sig11 <- y3_swept %*% Psi22ity3_swept + 1
  diag(Sig11) <- diag(Sig11) + 1

  Sig22iSig21 <- solve(Sig22, t(Sig12))
  y2mM2r <- (y2 - M2r)

  mu <- M1r + t(Sig22iSig21) %*% y2mM2r

  Sigcond <- Sig11 - Sig12 %*% Sig22iSig21

  Omegacond <- drop((Omega + t(y2mM2r) %*% solve(Sig22, y2mM2r))/(nu + Np))

  Sighat <- Sigcond * Omegacond

  Sighat <- (Sighat + t(Sighat))/2

  y1 <- mvnfast::rmvt(n = 1, mu = mu, sigma = Sighat, df = nu + Np)

  return(y1)

}

rcol_conditional_Mt <- function(y3y4 = NULL, nu = NULL, etap = NULL, etam = NULL, Psi22 = NULL, Psi22iPsi21 = NULL, Omega = NULL){

  ## This function draws a vector of discriminant observations if that particular discriminant is entirely missing

  y3y4_swept <- sweep(y3y4, 2, etap)

  M <- etam + y3y4_swept %*% Psi22iPsi21
  Sigcond <- y3y4_swept %*% solve(Psi22, t(y3y4_swept)) + 1
  diag(Sigcond) <- diag(Sigcond) + 1

  Sighat <- Sigcond*Omega/nu

  y1y2 <- mvnfast::rmvt(n = 1, mu = M, sigma = Sighat, df = nu)

  return(y1y2)

}

BayesECM_pred <- function(BayesECM_obj = NULL, Ytilde = NULL, verb = NULL,
                               thinning = NULL, eq_wts = NULL){


  la_p <- BayesECM_obj$MCMC$la_p
  a_p <- BayesECM_obj$MCMC$a_p
  Y <- BayesECM_obj$Y
  priors <- BayesECM_obj$priors


  N <- unname(sapply(Y, nrow))
  p <- max(sapply(Y, ncol))
  K <- length(Y)

  Kfull <- BayesECM_obj$data$Kfull
  if(length(Kfull) == 0){
    Kmissing <- 1:K
  }else{
    Kmissing <- (1:K)[-Kfull]
  }

  nu <- sapply(priors, function(X){X$nu})
  dof <- N + nu - p + 1

  Ntilde <- nrow(Ytilde)


  Ytilde_missing <- apply(Ytilde, 1, function(X){which(is.na(X))}, simplify = FALSE)

  umissing <- unique(Ytilde_missing)
  if(length(umissing) != 0){
  umissing <- umissing[sort(sapply(umissing, length), index.return = TRUE)$ix]
  }
  Ytilde_group <- rep(NA, times = Ntilde)

  for(i in 1:Ntilde){
    for(j in 1:length(umissing)){
      if(length(Ytilde_missing[[i]]) == length(umissing[[j]])){
        if(length(Ytilde_missing[[i]]) == 0){
          Ytilde_group[i] <- j
          break
          }else if(all(Ytilde_missing[[i]] == umissing[[j]])){
          Ytilde_group[i] <- j
          break
        }
      }
    }
  }



  group_ptilde <- p - sapply(umissing, length)

  ### MCMC

  predmean <- BayesECM_obj$MCMC$predmean
  SighatN <- BayesECM_obj$MCMC$SighatN
  Sighat <- matrix(NA, ncol = p, nrow = p)
  UT_template <- upper.tri(diag(p), diag = TRUE)
  LT_template <- lower.tri(diag(p))
  dens_template <- rep(NA, times = Ntilde)
  densmat <- matrix(0, ncol = K, nrow = Ntilde)

  dof <- unname(dof)
  Ytildemat <- unname(data.matrix(Ytilde))


  if(length(umissing) == 1){
    if(length(umissing[[1]]) == 0){
    #Ytilde has no missing data
    tpred <- tpred_nomissing

    }else{
      # All rows of Ytilde are missing the same discriminant
      tpred <- tpred_allmissing1

    }
  }else{
    if(length(umissing[[1]]) == 0){
      #Some rows of Ytilde are full, others are not
      tpred <- tpred_somemissing

    }else{
      # All rows of Ytilde are missing various discriminants
      # No rows of Ytilde are full

      tpred <- tpred_allmissing

    }
  }

  if(length(Kfull) != K){

    iters_use <- seq(from = 1, to = ncol(predmean[[Kmissing[1]]]), by = thinning)


  }else{
    iters_use <- NA
  }

  if(is.null(names(Y))){
    colnames(densmat) <- as.character(1:K)
  }else{
    colnames(densmat) <- names(Y)
  }




  if(eq_wts){
    pzfn <- equal_weights
  }else{
    pzfn <- unequal_weights
  }


  for(k in Kfull){

    Sighat[UT_template] <- SighatN[[k]]
    Sighat[LT_template] <- t(Sighat)[LT_template]
    Sighat <- Sighat/dof[k]

    muhat <- predmean[[k]]

    densmat[,k] <- tpred(Ytilde = Ytildemat, muhat = muhat, Sighat = Sighat, umissing = umissing,
                          Ytilde_group = Ytilde_group, dens_template = dens_template, dof = dof[k])


  }


  if(length(Kfull) == K){


    ldensmat <- log(densmat)

    pz <- pzfn(ldensmat = ldensmat, la_p = la_p)


  }else{

  for(i in iters_use){

    for(k in Kmissing){

      Sighat[UT_template] <- SighatN[[k]][,i]
      Sighat[LT_template] <- t(Sighat)[LT_template]
      Sighat <- Sighat/dof[k]

      muhat <- predmean[[k]][,i]

      densmat[,k] <- densmat[,k] + tpred(Ytilde = Ytildemat, muhat = muhat, Sighat = Sighat, umissing = umissing,
                            Ytilde_group = Ytilde_group, dens_template = dens_template, dof = dof[k])


    }


  }

    densmat[,Kmissing] <- densmat[,Kmissing]/length(iters_use)

    ldensmat <- log(densmat)

    pz <- pzfn(ldensmat = ldensmat, la_p = la_p)

    }

  pz <- as.data.frame(pz)
names(pz) <- names(BayesECM_obj$Y)

return(list(pz = pz, umissing = umissing, Ytilde_group = Ytilde_group, thinning = iters_use))

}

unequal_weights <- function(ldensmat = NULL, la_p = NULL){


  lwdens <- sweep(ldensmat, MARGIN = 2, la_p, "+")

  wdens <- exp(lwdens)
  lsumwdens <- log(rowSums(wdens))

  lpz <- lwdens - lsumwdens
  pz <- exp(lpz)

  return(pz)

}

equal_weights <- function(ldensmat = NULL, la_p = NULL){

  densmat <- exp(ldensmat)
  lsumdens <- log(rowSums(densmat))
  lpz <- ldensmat - lsumdens
  return(exp(lpz))

}

tpred_nomissing <- function(Ytilde = NULL, muhat = NULL, Sighat = NULL,
                            umissing = NULL, Ytilde_group = NULL, dens_template = NULL, dof = NULL){
  #no elements of Ytilde are missing

  dens_template <- mvnfast::dmvt(X = Ytilde, mu = muhat, sigma = Sighat, df = dof, log = FALSE)

  return(dens_template)

}

tpred_somemissing <- function(Ytilde = NULL, muhat = NULL, Sighat = NULL,
                              umissing = NULL, Ytilde_group = NULL, dens_template = NULL, dof = NULL){
  # Some rows of Ytilde have missing elements
  # Some rows of Ytilde are complete

  # dens_template is a previous allocation of memory for the return of this function


  dens_template[Ytilde_group == 1] <- mvnfast::dmvt(X = Ytilde[(Ytilde_group == 1), , drop = FALSE],
                                                   mu = muhat, sigma = Sighat, df =  dof, log = FALSE)


  for(j in 2:length(umissing)){


    dens_template[Ytilde_group == j] <- mvnfast::dmvt(X = Ytilde[(Ytilde_group == j), -umissing[[j]] , drop = FALSE],
                                                      mu = muhat[-umissing[[j]]],
                                                      sigma = Sighat[-umissing[[j]], -umissing[[j]]],
                                                      df =  dof, log = FALSE)


  }
  return(dens_template)
}

tpred_allmissing1 <- function(Ytilde = NULL, muhat = NULL, Sighat = NULL, dof = NULL,
                              umissing = NULL, Ytilde_group = NULL, dens_template = NULL){
  # All rows of Ytilde are missing the same elements

  dens_template <- mvnfast::dmvt(X = Ytilde[,-umissing[[1]] , drop = FALSE],
                                 mu = muhat[-umissing[[1]]], sigma = Sighat[-umissing[[1]], -umissing[[1]]],
                                 df = dof, log = FALSE)

  return(dens_template)

}

tpred_allmissing <- function(Ytilde = NULL, muhat = NULL, Sighat = NULL, dof = NULL,
                             umissing = NULL, Ytilde_group = NULL, dens_template = NULL){
  # All rows of Ytilde have a missing element, but what is missing varries

  for(j in 1:length(umissing)){

    dens_template[Ytilde_group == j] <- mvnfast::dmvt(X =  Ytilde[(Ytilde_group == j), -umissing[[j]] , drop = FALSE],
                                                      mu = muhat[-umissing[[j]]],
                                                      sigma = Sighat[-umissing[[j]], -umissing[[j]]],
                                                      df = dof, log = FALSE)

  }
  return(dens_template)

}


#' Summary of Unlabeled Event Categorization
#'
#' Tabulates results from the [predict.BayesECM()] function for quick analysis.
#'
#' @param object an object of `class` `"BayesECMpred"` obtained as the output from the [predict.BayesECM()] function.
#' @param index integer stipulating the event of interest.  Value corresponds to the row index of `Ytilde` previously supplied to [predict.BayesECM()]
#' @param category integer for the index of the category of interest for hypothesis testing.  Alternatively, a character string naming the category of interest can be provided.
#' @param C square matrix of dimension 2, providing loss values to be used in hypothesis testing.  See Details.
#' @param ... not used
#'
#' @return Prints a summary including probability of each category for the event stipulated by `index`, minimum expected loss for binary categorization, and probability of a-typicality of the event for the category specified by `category`.
#' @export
#'
#' @details
#'
#' ## Expected Loss
#'
#' [summary.BayesECMpred()] prints expected loss for the binary hypothesis stipulated by `category`.  Expected loss is calculated using the loss matrix specified with argument `C`.  The default values for `C` result in 0-1 loss being used.  Format details for the loss matrix can be found in `vignette("syn-data-code")`.
#'
#' ## Typicality
#'
#' Typicality indices are used in [cECM_decision()] as part of the decision criteria.  Here, we have adapted typicality indices for use with a Bayesian ECM model for outlier detection, when a new observation may not be related to the categories used for training.  Probability of the p-value being less than a significance level of 0.05 is reported.  If no missing data is used for training, this probability is either 0 or 1.
#'
#' @examples
#'
#' csv_use <- "good_training.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' training_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' trained_model <- BayesECM(Y = training_data, BT = c(10,1000))
#'
#' csv_use <- "good_newdata.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' new_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' bayespred <- predict(trained_model,  Ytilde = new_data)
#'
#'
#' @method summary BayesECMpred
#'
summary.BayesECMpred <- function(object, index = 1, category = 1, C = 1 - diag(2), ...){

  X <- as.data.frame(object$epz[index , , drop = FALSE])
  cat_names <- colnames(X)

  if(is.numeric(category)){
    if(category > length(cat_names)){
      stop("Integer value for the argument 'category' provided.  The value provided selects the category index and must be less than or equal to the number of total event categories.")
    }
    cat_index <- category
    category <- cat_names[category]

  }else if(is.character(category)){
    cat_index <- which(cat_names == category)
    if(length(cat_index) == 0){
      stop("Supplied argument 'category' is not found within the categories used to generate 'object'.  Check for correct spelling of the 'category' argument.")
    }
  }else{
    stop("Supplied argument 'category' must be an integer or character string which provides the index or names the category of interest.  Typically 'category' corresponds to the name used to specify detonations.")
  }


  cat_samples <- X[category]

  else_samples <- X[,-which(category == dimnames(X)[[2]]) ,drop = FALSE]

  else_samples <- rowSums(else_samples)

  Xall <- cbind(X,else_samples)

  reorder <- c(cat_index, ncol(Xall), (1:(ncol(Xall)-1))[-cat_index])
  Xall <- Xall[,reorder]

  # First data.frame

  stat.sum <- data.frame(matrix(NA, ncol = 1, nrow = ncol(Xall)))
  names(stat.sum) <- c("E")

  # Expectation

  stat.sum$E <- unlist(unname(Xall))

  rownames(stat.sum) <- c(category, paste0("Not ", category), colnames(X)[-cat_index])
  stat.sum$E <- round(stat.sum$E, digits = 3)
  ellps <- data.frame(E = rep(".....", times = 1))
  rownames(ellps) <- "...Itemized Categories..."
  stat.sum <- rbind(stat.sum[1:2, , drop = FALSE], ellps , stat.sum[3:nrow(stat.sum), , drop = FALSE])


  # Print stat.sum data.frame

  cat("\n\n")
  cat("Summary Statistics")
  print(knitr::kable(stat.sum, format = "rst", digits = 2, row.names = TRUE, col.names = c("Expected Probability of Category")))

  # New data.frame

  loss_df <- data.frame(matrix(NA, ncol = 2, nrow = 2))
  names(loss_df) <- c("EL", "MinEL")
  # Expected Loss
  # Arrow Pointing to minimum expected loss

  Ep <- as.numeric(stat.sum[1:2, 1])
  a1 <- sum(C[,1]*Ep)
  a2 <- sum(C[,2]*Ep)

  loss_df$EL <- c(a1,a2)
  loss_df$MinEL[which.min(loss_df$EL)] <- "<-"
  loss_df$MinEL[is.na(loss_df$MinEL)] <- " "

  rownames(loss_df) <- rownames(stat.sum)[1:2]

  percent_crit_nuem <- (C[2,2] - C[2,1])
  percent_crit_denom <- C[1,1] - C[2,1] + C[2,2] - C[1,2]

  if(percent_crit_denom < 0){
    equality_sign <- ">"
  }else{
    equality_sign <- "<"
  }

  cat("\n\n")
cat(paste0("Decision Criterion: select ", category, " if E[p(", category, ")] ", equality_sign, " ", round(percent_crit_nuem/percent_crit_denom, digits = 2), ",\nequivalent to minimum expected loss"))
  print(knitr::kable(loss_df, format = "rst", row.names = TRUE, col.names = c("Expected Loss", "Minimum Expected Loss")))

  ## typicality index

  Kfull <- object$BayesECMfit$data$Kfull
  K <- length(object$BayesECMfit$Y)

  if(length(Kfull) == 0){
    Kmissing <- 1:K
    iters <- ncol(object$BayesECMfit$MCMC$Ybar[[Kmissing[1]]])
  }else if(all(Kfull %in% (1:K))){
    Kmissing <- integer(0)
    iters <- NA
  }else{
    Kmissing <- (1:K)[which(!((1:K) %in% Kfull))]
    iters <- ncol(object$BayesECMfit$MCMC$Ybar[[Kmissing[1]]])
  }

  thinning <- object$thinning

  Ytilde <- unlist(unname(object$Ytilde[index,]))
  Ytilde_missing <- which(is.na(Ytilde))
  Ytilde <- Ytilde[-Ytilde_missing]
  ptilde <- length(Ytilde)
  pvals <- matrix(NA, ncol = K, nrow = 1)

  N <- sapply(object$BayesECMfit$Y, nrow)
  nu <- sapply(object$BayesECMfit$priors, function(X){
    X$nu
  })
  p <- ncol(object$BayesECMfit$Y[[1]])
  Sighat <- matrix(NA, ncol = p, nrow = p)
  UT_template <- upper.tri(diag(p), diag = TRUE)
  LT_template <- lower.tri(diag(p))
  Sighat <- matrix(NA, ncol = p, nrow = p)
  dof <- N + nu - p + 1

  typicality_probs <- rep(NA, times = K)
  alphatilde <- 0.05

  if(length(Kfull) != 0){

    for(k in Kfull){

      Sighat[UT_template] <- object$BayesECMfit$MCMC$SighatN[[k]]
      Sighat[LT_template] <- t(Sighat)[LT_template]

      muhat <- object$BayesECMfit$MCMC$predmean[[k]]

      if(ptilde == p){
        muuse <- muhat
        Siguse <- Sighat
      }else{
        muuse <- muhat[-Ytilde_missing, drop = FALSE]
        Siguse <- Sighat[-Ytilde_missing,-Ytilde_missing, drop = FALSE]
      }

      y <- Ytilde - muuse
      Siy <- solve(Siguse,y)
      tySiy <- t(y) %*% Siy

      pval <- stats::pf(tySiy*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)

      typicality_probs[k] <- as.numeric(pval < alphatilde)


    }

  }


  if(length(Kmissing) != 0){

    for(k in Kmissing){
      pval_samples <- rep(NA, times = iters)
      for(i in thinning){

        Sighat[UT_template] <- object$BayesECMfit$MCMC$SighatN[[k]][,i]
        Sighat[LT_template] <- t(Sighat)[LT_template]

        muhat <- object$BayesECMfit$MCMC$predmean[[k]][,i]

        if(ptilde == p){
          muuse <- muhat
          Siguse <- Sighat
        }else{
          muuse <- muhat[-Ytilde_missing, drop = FALSE]
          Siguse <- Sighat[-Ytilde_missing,-Ytilde_missing, drop = FALSE]
        }

        y <- Ytilde - muuse
        Siy <- solve(Siguse,y)
        tySiy <- t(y) %*% Siy

        pval <- stats::pf(tySiy*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)

        pval_samples[i] <- pval < alphatilde

      }


      typicality_probs[k] <- mean(pval_samples, na.rm = TRUE)

    }

  }

  cat("\n\n")
  cat(paste0("Probaility of a-typicality of category ", category, " is ", round(typicality_probs[cat_index], digits = 3), " for a significance level of 0.05"))

 }


#' B-ECM performance metrics
#'
#' Outputs batch performance metrics of decisions using the output of [predict.BayesECM()] when the true category of testing events are known.  Can be used for empirically comparing different model fits.
#'
#' @param bayes_pred An object of class `"BayesECMpred"` returned from the [predict.BayesECM()] function.  Data and order of data provided to [predict.BayesECM()] must match the order of the `cat_truth` argument.
#' @param vic Character string, indicating the "Very Important Category" (`vic`) used for calculating categorization metrics.  The return of `becm_decision()` provides information on if new observations should be categorized into `vic` or the remainder of the training categories grouped together.
#' @param cat_truth Vector of the same length as `nrow(bayes_pred$Ytilde)`, where `Ytilde` was previously supplied to the [predict.BayesECM()] function.  Used for accuracy calculations of a `BayesECM` fit and decision criteria when the true category of the new observations are known.
#' @param alphatilde Numeric scalar between 0 and 1 used for the significance level of typicality indices.
#' @param pn Logical.  Acronym for "false Positive, false Negative". `pn == FALSE` indicates that only accuracy of categorizations should be returned, while `pn == TRUE` indicates that false positives and false negatives should be returned in addition to accuracy.
#' @param C Square matrix of dimension 2 or the number of categories.  Used as the loss function in the decision theoretic framework.  The default is 0-1 loss for binary categorization.
#' @param rej `data.frame` of rejection via typicality index retrieved from a previous call to `becm_decision()`.  Useful for saving computation time when comparing the results from different supplied matrices for argument `C`, while using a constant value of `alphatilde`.
#'
#' @return A list of two data frames of logicals.  The rows in each data frame correspond to the rows in `bayes_pred$Ytilde`.  The first data frame, named `results`, has three columns named `correct`, `fn`, and `fp`. The `results` column indicates if the categorization is correct.  `fn` and `fp` stand for false negatives and false positives respectively.  `fn` and `fp` are found under binary categorization.  Values of `NA` are returned when a false positive or false negative is not relevant.  The second data frame, named `rej`, indicates the rejection of each new observation from each category via typicality index.  One can use `rej` to inspect if the typicality index played a role in categorization, or to supply to another call of `becm_decision`.
#'
#' @details The matrix `C` specifies the loss for a set of categorization actions for a single new observation \eqn{\tilde{y}_{\tilde{p}}} given the probability of \eqn{\tilde{y}_{\tilde{p}}} belonging to each of \eqn{K} categories.  Actions are specified as the columns, and the event category random variables are specified as the rows.  See the vignette \code{vignette("syn-data-code", package = "ezECM")} for more mathematical details.
#'
#' The dimension of matrix `C` specifies the categorization type.  A dimension of 2 is binary categorization, with the first column and row always corresponding to the category chosen as the `vic` argument.  Otherwise, when the dimension of `C` is equal to the number of categories, the indices of the rows and columns of `C` are in the same order as the categories listed for `names(bayes_pred$BayesECMfit$Y)`.
#'
#' @examples
#'
#' csv_use <- "good_training.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' training_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' trained_model <- BayesECM(Y = training_data, BT = c(10,1000))
#'
#' csv_use <- "good_newdata.csv"
#' file_path <- system.file("extdata", csv_use, package = "ezECM")
#' new_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' bayespred <- predict(trained_model,  Ytilde = new_data)
#'
#' accuracy <- becm_decision(bayes_pred = bayespred, vic = "explosion", cat_truth = new_data$event)
#'
#'
#' @export
#'
becm_decision <- function(bayes_pred = NULL, vic = NULL, cat_truth = NULL, alphatilde = 0.05, pn = TRUE, C = matrix(c(0,1,1,0), ncol = 2), rej = NULL){

  if(!inherits(bayes_pred, "BayesECMpred")){
    stop("Argument 'bayes_pred' must be the output of the predict.BayesECM function.")
  }

  if(length(alphatilde) != 1 | !is.numeric(alphatilde)){
    stop("Argument 'alphatilde' must be a numeric of length 1.")
  }

  if(alphatilde > 1 | alphatilde < 0){
    stop("Argument 'alphatilde' must be a value betweeen 0 and 1.")
  }

  if(length(vic) != 1 ){
    stop("Argument 'vic' must be a vector of length 1.")
  }

  if(!is.character(vic)){
    warning("Attempting to coerce non character value supplied to the 'vic' argument to a character.", immediate. = TRUE)
    vic <- as.character(vic)
  }

  if(!(vic %in% names(bayes_pred$BayesECMfit$Y))){
    stop("Argument 'vic' must be one of the category names stipulated in training the BayesECM model.")
  }

  if(length(cat_truth) != nrow(bayes_pred$Ytilde)){
    stop("Argument 'cat_truth' must have the same length as nrow(bayes_pred$Ytilde).")
  }

  if(!all(cat_truth %in% names(bayes_pred$BayesECMfit$Y))){
    stop("Values of argument 'cat_truth' must all be present in the original categories used for training.")
  }

  if(!(nrow(C) %in% c(2, length(names(bayes_pred$BayesECMfit$Y)))) | !is.matrix(C) | !is.numeric(C) | nrow(C) != ncol(C)){
    stop("Argument 'C' must be a numeric square matrix.  The dimension of C must be either equal to 2 or the number of categories used in the training data set.")
  }


  if(length(bayes_pred$BayesECMfit$data$Kfull) == length(bayes_pred$BayesECMfit$Y)){
    out <- becm_decision_fulltraining(bayes_pred = bayes_pred, alphatilde = alphatilde, vic = vic, cat_truth = cat_truth, pn = pn, C = C, rej = rej)
  }else{

    out <- becm_decision_missingtraining(bayes_pred = bayes_pred, alphatilde = alphatilde, vic = vic, cat_truth = cat_truth, pn = pn, C = C, rej = rej)

  }

  return(out)


}


becm_decision_fulltraining <- function(bayes_pred = NULL, alphatilde = NULL, vic = NULL, cat_truth = NULL, pn = NULL, C = NULL, rej = NULL){

  K <- length(bayes_pred$BayesECMfit$Y)

  Ytilde <- bayes_pred$Ytilde

  pvals <- matrix(NA, ncol = K, nrow = nrow(bayes_pred$Ytilde))
  N <- sapply(bayes_pred$BayesECMfit$Y, nrow)
  nu <- sapply(bayes_pred$BayesECMfit$priors, function(X){
    X$nu
  })
  p <- ncol(bayes_pred$BayesECMfit$Y[[1]])
  Sighat <- matrix(NA, ncol = p, nrow = p)
  UT_template <- upper.tri(diag(p), diag = TRUE)
  LT_template <- lower.tri(diag(p))
  Sighat <- matrix(NA, ncol = p, nrow = p)
  dof <- N + nu - p + 1


  if(nrow(C) == 2){

  pAB <- t(apply(bayes_pred$epz, 1, function(X,v){
      x <- rep(NA, times = 2)
      x[1] <- X[which(names(X) == v)]
      x[2] <- sum(X[which(names(X) != v)])
      return(x)
  }, v = vic, simplify = TRUE))


  Cdenom <- C[1,2] - C[2,2] + C[2,1] - C[1,1]
  Cratio <- (C[2,1] - C[2,2])/Cdenom

  pcat <- rep("b", times = nrow(pAB))
  if(Cdenom > 0){
    pcat[pAB[,1] > Cratio] <- "a"
  }else{
    pcat[pAB[,1] < Cratio] <- "a"
  }
  }else{


    Zdist_missing <- bayes_pred$epz
    Eloss <-  as.matrix(Zdist_missing) %*% C
    minEloss <- apply(Eloss, 1, which.min)
    vicindex <- which(names(Zdist_missing[1,]) == vic)
    vicminEloss <- (minEloss == which(names(Zdist_missing[1,]) == vic))
    pcat <- rep("b", times = nrow(Ytilde))
    pcat[vicminEloss] <- "a"
    minEloss_char <- names(Zdist_missing[1,])[minEloss]

}

  if(is.null(rej)){
    rej <- matrix(NA, ncol = K, nrow = nrow(Ytilde))
  if(length(bayes_pred$umissing) != 1){
    for(k in 1:K){

      pvals <- rep(NA, times = nrow(Ytilde))

      Sighat[UT_template] <- bayes_pred$BayesECMfit$MCMC$SighatN[[k]]
      Sighat[LT_template] <- t(Sighat)[LT_template]

      muhat <- bayes_pred$BayesECMfit$MCMC$predmean[[k]]

      for(i in sort(unique(bayes_pred$Ytilde_group))){

        if(length(bayes_pred$umissing[[1]]) == 0 & i == 1){

          ind <- which(bayes_pred$Ytilde_group == 1)

          ptilde <- p - length(bayes_pred$umissing[[i]])

          pvals[ind] <- apply(Ytilde[ind, , drop = FALSE], 1, function(X,S, ptil, m){
            y <- (X - m)
            Siy <- solve(S,y)
            tySiy <- drop(t(y) %*% Siy)
            return(tySiy/ptil)
          }, S = Sighat, ptil = ptilde, m = muhat)

          pvals[ind] <- stats::pf(pvals[ind]*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)

        }else{
          ind <- which(bayes_pred$Ytilde_group == i)

          ptilde <- p - length(bayes_pred$umissing[[i]])

          muuse <- muhat[-bayes_pred$umissing[[i]], drop = FALSE]
          Siguse <- Sighat[-bayes_pred$umissing[[i]], -bayes_pred$umissing[[i]], drop = FALSE]

          pvals[ind] <- apply(Ytilde[ind,-bayes_pred$umissing[[i]], drop = FALSE], 1, function(X,S, ptil, m){
            y <- (X - m)
            Siy <- solve(S,y)
            tySiy <- t(y) %*% Siy
            #tySiy <- (tySiy + t(tySiy))/2
            return(tySiy/ptil)
          }, S = Siguse, ptil = ptilde, m = muuse)

          pvals[ind] <- stats::pf(pvals[ind]*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)
        }
      }
      rej[,k] <- pvals < alphatilde
    }
  }else{
    if(length(bayes_pred$umissing[[1]]) == 0){

      for(k in 1:K){

        Sighat[UT_template] <- bayes_pred$BayesECMfit$MCMC$SighatN[[k]]
        Sighat[LT_template] <- t(Sighat)[LT_template]

        muhat <- bayes_pred$BayesECMfit$MCMC$predmean[[k]]

        pvals <- apply(Ytilde, 1, function(X,S, ptilde, m){
          y <- (X - m)
          Siy <- solve(S,y)
          tySiy <- t(y) %*% Siy
          return(tySiy/ptilde)
        }, S = Sighat, ptilde = p, m = muhat)


        pvals <- stats::pf(pvals*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)


        rej[,k] <- pvals < alphatilde

      }


    }else{

      ptilde <- p - length(bayes_pred$umissing[[1]])

      for(k in 1:K){

      Sighat[UT_template] <- bayes_pred$BayesECMfit$MCMC$SighatN[[k]]
      Sighat[LT_template] <- t(Sighat)[LT_template]

      muhat <- bayes_pred$BayesECMfit$MCMC$predmean[[k]]

      muuse <- muhat[-bayes_pred$umissing[[1]], drop = FALSE]
      Siguse <- Sighat[-bayes_pred$umissing[[1]], -bayes_pred$umissing[[1]], drop = FALSE]

      pvals <- apply(Ytilde[,-bayes_pred$umissing[[1]], drop = FALSE], 1, function(X,S, ptil, m){
        y <- (X - m)
        Siy <- solve(S,y)
        tySiy <- t(y) %*% Siy
        return(tySiy/ptil)
      }, S = Siguse, ptil = ptilde, m = muuse)

      pvals <- stats::pf(pvals*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)


      rej[,k] <- pvals < alphatilde

      }
      }
  }

}


  for(i in 1:length(pcat)){
    if(pcat[i] == "a"){
      if(rej[i,which(colnames(bayes_pred$epz) == vic)]){
        pcat[i] <- "b"
        if(nrow(C) > 2){
          minEloss_char[i] <- "outlier"
        }
      }
    }
  }

  true_group <- cat_truth
  true_group[which(cat_truth != vic)] <- "b"
  true_group[which(cat_truth == vic)] <- "a"

  if(pn){

    out <- data.frame(matrix(NA, ncol = 3, nrow = nrow(Ytilde)))
    names(out) <- c("correct", "fn", "fp")

    if(nrow(C) == 2){
      out$correct <- pcat == true_group
      out$fn[cat_truth == vic] <- !(out$correct[cat_truth == vic])
      out$fp[cat_truth != vic] <- !(out$correct[cat_truth != vic])
    }else{

      out$correct <- minEloss_char == cat_truth
      bicat <-  pcat == true_group
      out$fn[cat_truth == vic] <- !(bicat[cat_truth == vic])
      out$fp[cat_truth != vic] <- !(bicat[cat_truth != vic])

    }



  }else{

    if(nrow(C) == 2){
      out <- pcat == true_group

    }else{

      out <- minEloss_char == cat_truth

    }
  }

  rej <- data.frame(rej)
  names(rej) <- names(bayes_pred$BayesECMfit$Y)

  return(list(results = out, rej = rej))

}

becm_decision_missingtraining <- function(bayes_pred = NULL, alphatilde = NULL, vic = NULL, cat_truth = NULL, pn = NULL, C = NULL, rej = NULL){

  K <- length(bayes_pred$BayesECMfit$Y)

  Kfull <- bayes_pred$BayesECMfit$data$Kfull
  if(length(Kfull) == 0){
    Kmissing <- 1:K
  }else{
    Kmissing <- (1:K)[-Kfull]
  }

  iters <- ncol(bayes_pred$BayesECMfit$MCMC$SighatN[[Kmissing[1]]])
  iters_use <- bayes_pred$thinning

  Ytilde <- bayes_pred$Ytilde


  N <- sapply(bayes_pred$BayesECMfit$Y, nrow)
  nu <- sapply(bayes_pred$BayesECMfit$priors, function(X){
    X$nu
  })
  p <- ncol(bayes_pred$BayesECMfit$Y[[1]])
  Sighat <- matrix(NA, ncol = p, nrow = p)
  UT_template <- upper.tri(diag(p), diag = TRUE)
  LT_template <- lower.tri(diag(p))
  Sighat <- matrix(NA, ncol = p, nrow = p)
  dof <- N + nu - p + 1

  Zdist_missing <-bayes_pred$epz

  if(nrow(C) == 2){
    pAB <- t(apply(Zdist_missing, 1, function(X, v){
      x <- rep(NA, times = 2)
      x[1] <- X[which(names(X) == v)]
      x[2] <- sum(X[which(names(X) != v)])
      return(x)
    }, v = vic))


    Cdenom <- C[1,2] - C[2,2] + C[2,1] - C[1,1]
    Cratio <- (C[2,1] - C[2,2])/Cdenom
    pcat <- rep("b", times = nrow(pAB))
    if(Cdenom > 0){
      pcat[pAB[,1] > Cratio] <- "a"
    }else{
      pcat[pAB[,1] < Cratio] <- "a"
    }
  }else{

    Eloss <- as.matrix(Zdist_missing) %*% C
    minEloss <- apply(Eloss, 1, which.min)

    vicindex <- which(names(Zdist_missing[1,]) == vic)

    vicminEloss <- (minEloss == which(names(Zdist_missing[1,]) == vic))

    pcat <- rep("b", times = nrow(Ytilde))
    pcat[vicminEloss] <- "a"

    minEloss_char <- names(Zdist_missing[1,])[minEloss]

  }


  if(is.null(rej)){
    rej <- matrix(NA, ncol = K, nrow = nrow(Ytilde))

    for(k in Kfull){

      pvals <- rep(NA, times = nrow(Ytilde))

      Sighat[UT_template] <- bayes_pred$BayesECMfit$MCMC$SighatN[[k]]
      Sighat[LT_template] <- t(Sighat)[LT_template]

      muhat <- bayes_pred$BayesECMfit$MCMC$predmean[[k]]

      for(j in sort(unique(bayes_pred$Ytilde_group))){

        if(length(bayes_pred$umissing[[1]]) == 0 & j == 1){

          ind <- which(bayes_pred$Ytilde_group == 1)

          ptilde <- p - length(bayes_pred$umissing[[j]])

          pvals[ind] <- apply(Ytilde[ind, , drop = FALSE], 1, function(X,S, ptil, m){
            y <- (X - m)
            Siy <- solve(S,y)
            tySiy <- drop(t(y) %*% Siy)
            return(tySiy/ptil)
          }, S = Sighat, ptil = ptilde, m = muhat)

          pvals[ind] <- stats::pf(pvals[ind]*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)

        }else{
          ind <- which(bayes_pred$Ytilde_group == j)

          ptilde <- p - length(bayes_pred$umissing[[j]])

          muuse <- muhat[-bayes_pred$umissing[[j]], drop = FALSE]
          Siguse <- Sighat[-bayes_pred$umissing[[j]], -bayes_pred$umissing[[j]], drop = FALSE]

          pvals[ind] <- apply(Ytilde[ind,-bayes_pred$umissing[[j]], drop = FALSE], 1, function(X,S, ptil, m){
            y <- (X - m)
            Siy <- solve(S,y)
            tySiy <- t(y) %*% Siy
            return(tySiy/ptil)
          }, S = Siguse, ptil = ptilde, m = muuse)

          pvals[ind] <-stats::pf(pvals[ind]*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)
        }
      }
      rej[,k] <- pvals < alphatilde
    }
    for(k in Kmissing){

      pvals <- matrix(NA, nrow = nrow(Ytilde), ncol = iters_use[length(iters_use)])

      for(i in iters_use){

        Sighat[UT_template] <- bayes_pred$BayesECMfit$MCMC$SighatN[[k]][,i]
        Sighat[LT_template] <- t(Sighat)[LT_template]

        muhat <- bayes_pred$BayesECMfit$MCMC$predmean[[k]][,i]

        for(j in sort(unique(bayes_pred$Ytilde_group))){

          if(length(bayes_pred$umissing[[1]]) == 0 & j == 1){

            ind <- which(bayes_pred$Ytilde_group == 1)

            ptilde <- p - length(bayes_pred$umissing[[j]])

            pvals[ind,i] <- apply(Ytilde[ind, , drop = FALSE], 1, function(X,S, ptil, m){
              y <- (X - m)
              Siy <- solve(S,y)
              tySiy <- drop(t(y) %*% Siy)
              return(tySiy/ptil)
            }, S = Sighat, ptil = ptilde, m = muhat)

            pvals[ind,i] <- stats::pf(pvals[ind,i]*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)

          }else{
            ind <- which(bayes_pred$Ytilde_group == j)

            ptilde <- p - length(bayes_pred$umissing[[j]])

            muuse <- muhat[-bayes_pred$umissing[[j]], drop = FALSE]
            Siguse <- Sighat[-bayes_pred$umissing[[j]], -bayes_pred$umissing[[j]], drop = FALSE]

            pvals[ind,i] <- apply(Ytilde[ind,-bayes_pred$umissing[[j]], drop = FALSE], 1, function(X,S, ptil, m){
              y <- (X - m)
              Siy <- solve(S,y)
              tySiy <- t(y) %*% Siy
              return(tySiy/ptil)
            }, S = Siguse, ptil = ptilde, m = muuse)

            pvals[ind,i] <- stats::pf(pvals[ind,i]*dof[k], df1 = ptilde, df2 = dof[k], log.p = FALSE, lower.tail = FALSE)
          }
        }
      }



      pvals <- apply(pvals, 1, function(X,a){
        X <- X[!is.na(X)]
        FX <- sum(X < a)/length(X)
        return(FX)
      }, a = alphatilde)


      rej[,k] <- pvals > 0.5
    }
  }


  for(i in 1:length(pcat)){
    if(pcat[i] == "a"){
      if(rej[i,which(names(bayes_pred$BayesECMfit$Y) == vic)]){
        pcat[i] <- "b"
        if(nrow(C) > 2){
          minEloss_char[i] <- "outlier"
        }
      }
    }
  }

  true_group <- cat_truth
  true_group[which(cat_truth != vic)] <- "b"
  true_group[which(cat_truth == vic)] <- "a"

  if(pn){

    out <- data.frame(matrix(NA, ncol = 3, nrow = nrow(Ytilde)))
    names(out) <- c("correct", "fn", "fp")

    if(nrow(C) == 2){
      out$correct <- pcat == true_group
      out$fn[cat_truth == vic] <- !(out$correct[cat_truth == vic])
      out$fp[cat_truth != vic] <- !(out$correct[cat_truth != vic])
    }else{

      out$correct <- minEloss_char == cat_truth
      bicat <-  pcat == true_group
      out$fn[cat_truth == vic] <- !(bicat[cat_truth == vic])
      out$fp[cat_truth != vic] <- !(bicat[cat_truth != vic])

    }



  }else{

    if(nrow(C) == 2){
      out <- pcat == true_group

    }else{

      out <- minEloss_char == cat_truth

    }
  }

  rej <- data.frame(rej)
  names(rej) <- names(bayes_pred$BayesECMfit$Y)
  return(list(results= out, rej = rej))

}
