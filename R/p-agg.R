#' Multiple Discriminant Analysis
#'
#' Fits a regularized discriminant analysis model to labeled training data and generates an aggregate p-value for categorizing newly obtained data.
#'
#' Details on regularized discriminant analysis (RDA) can be found in \insertCite{friedman1989regularized;textual}{ezECM}.  Details on related implementation found in \insertCite{anderson2007mathematical;textual}{ezECM}.
#'
#' @param x Either a `data.frame` of training data or the `list` output from previous use of the `cECM` function.  If a `data.frame` is supplied, each row contains one point of training data.  Additionally, one column, with the name of `event`, must contain labels for the event of each row.  The remaining columns contain the observed `ncol(x) - 1` discriminant related p-value data.
#' @param newdata a `data.frame` containing new data for event categorization.  Must use a subset of the discriminants used in the training data.  If for a particular event a certain discriminant is unavailable, specifying `NA` in place of the p-value will allow for calculation of the aggregate p-value.
#' @param rda_params a `list` of arguments passed to the [klaR::rda()] function.  If arguments `rda_params$x` or `rda_params$grouping` are supplied, they are ignored.
#' @param transform Logical indicating if the supplied p-values should be transformed by the function \eqn{2/\pi \times \mathrm{asin}\sqrt{X}}.  Ignored if a `list` is supplied as the argument `x`.
#'
#' @return A list.  Any returned objects contain a list element indicating the value of `transform` supplied to the `cECM` function call, as well as a [klaR::rda()] object related to relevant training data.  In addition if `newdata` argument is supplied, the returned list contains a `data.frame` specifying aggregate p-values for each new event (rows) for related event category (columns).
#'
#' @export
#'
#' @importFrom Rdpack reprompt
#' @import stats
#' @import klaR
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#'
#' x <- pval_gen(sims = 20, pwave.arrival = list(optim.starts = 5))
#' s <- sample(1:20, size = 2)
#'
#' newdata <- x[s,]
#' newdata <- newdata[,-which(names(newdata) == "event")]
#'
#' x <- x[-s,]
#'
#' pval_cat <- cECM(x = x, transform = TRUE)
#'
#' pval_cat <- cECM(x = pval_cat, newdata = newdata)
#'
#'
cECM <- function(x, newdata = NULL, rda_params = NULL, transform = TRUE){



  ## TO DO
  ### update documentation regarding function output
  ###  Ensure newdata columns are in the same order as the training data and that the names match

  eps <- sqrt(.Machine$double.eps)

  if(!is.null(rda_params)){
    rda_params$x <- NULL
    rda_params$grouping <- NULL
  }


  if(inherits(x, "data.frame") ){

    if(!any(names(x) == "event")){
      warning("A column of x naming the event class for each row must be stipulated.", immediate. = TRUE)
    }

    if(transform == TRUE){
      x.data <- which(!(names(x) == "event"))
      x[, x.data] <- apply(x[, x.data], 2, function(X){
        return((2/pi)*asin(sqrt(X)))
        #return(log(X) -log(1-X))
      })
    }

    rda_fit <- klaR::rda(x = x[,-which(names(x) == "event"), drop = FALSE], grouping = x$event)#do.call(klaR::rda, c(list(x = x[,-which(names(x) == "event")], grouping = x$event), rda_params))


    if(!is.null(newdata)){

      if(transform == TRUE){
        newdata[, which(names(newdata) != "event")] <- apply(newdata[, which(names(newdata) != "event")], 2, function(X){
          return((2/pi)*asin(sqrt(X)))
          #return(log(X) -log(1-X))
        })
      }

      newdata.save <- newdata
      newdata$event <- NULL

      if(!all(names(newdata) %in% names(x[,-which(names(x) == "event")]))){
        stop("Discriminant names for newdata are inconsistent with discriminant names of the training data.  This must be corrected before using the cECM() function.")
      }

      if(any(names(newdata) != names(x[,-which(names(x) == "event")]))){
        m <- match(names(newdata), names(x[,-which(names(x) == "event")]))
        newdata <- newdata[,m]
        newdata.save <- newdata.save[, c(m, ncol(newdata.save))]
      }



    p <- length(names(x))-1
    K <- length(unique(x$event))

    cECM_out <- data.frame(matrix(NA, ncol = K, nrow = nrow(newdata)))
    names(cECM_out) <- rda_fit$classes#names(x)[-which(names(x) == "event")]
    llik <- cECM_out

    for(k in 1:K){

      Sig <- rda_fit$covariances[, , k]
      Sig <- Sig + diag(nrow(Sig))*eps

      mu <- rda_fit$means[,k]

      cECM_vec <- llik_vec <- rep(NA, times = nrow(newdata))

      for(i in 1:nrow(newdata)){

        Y <- unlist(unname(newdata[i,]))

        na.index <- which(!is.na(Y))
        Y <- Y[na.index]
        mu.use <- mu[na.index]

        Sig.use <- Sig[na.index, na.index, drop = FALSE]



        ldet <- determinant(Sig.use, logarithm = TRUE)$modulus

        chi2.stat <- drop(crossprod(Y - mu.use, solve(Sig.use,(Y - mu.use))))



        cECM_vec[i] <- 1-stats::pchisq(chi2.stat, df = length(Y))
        llik_vec[i] <- chi2.stat + ldet

      }

      cECM_out[,k] <- cECM_vec
      llik[,k] <- llik_vec

    }

    }




  }else if(inherits(x, "cECM") ){

    rda_fit <- x$rda_fit

    if(!is.null(newdata)){


      if(x$transform == TRUE){
        newdata[, which(names(newdata) != "event")] <- apply(newdata[, which(names(newdata) != "event")], 2, function(X){
          return((2/pi)*asin(sqrt(X)))
        })
      }

      newdata.save <- newdata
      newdata$event <- NULL

      if(!all(names(newdata) %in% names(x$x[,-which(names(x$x) == "event")]))){
        stop("Discriminant names for newdata are inconsistent with discriminant names of the training data.  This must be corrected before using the cECM() function.")
      }

      if(all(names(newdata) != names(x$x[,-which(names(x$x) == "event")]))){
        m <- match(names(newdata), names(x$x[,-which(names(x$x) == "event")]))
        newdata <- newdata[,m]
        newdata.save <- newdata.save[, c(m, ncol(newdata.save))]
      }

      K <- length(rda_fit$classes)

      p <- ncol(x$x) - 1

      cECM_out <- data.frame(matrix(NA, ncol = K, nrow = nrow(newdata)))
      names(cECM_out) <- rda_fit$classes
      llik <- cECM_out

      for(k in 1:K){

        Sig <- rda_fit$covariances[, , k]
        Sig <- Sig + diag(nrow(Sig))*eps



        mu <- rda_fit$means[,k]

        cECM_vec <- llik_vec <- rep(NA, times = nrow(newdata))

        for(i in 1:nrow(newdata)){

          Y <- unlist(unname(newdata[i,]))

          na.index <- which(!is.na(Y))
          Y <- Y[na.index]
          mu.use <- mu[na.index]

          A <- diag(length(mu))
          A <- A[na.index, , drop = FALSE]

          Sig.use <- A %*% Sig %*% t(A)

          ldet <- determinant(Sig.use, logarithm = TRUE)$modulus

          chi2.stat <- drop(crossprod(Y - mu.use, solve(Sig.use,(Y - mu.use))))



          cECM_vec[i] <- 1-stats::pchisq(chi2.stat, df = length(Y))
          llik_vec[i] <- chi2.stat + ldet


        }

        cECM_out[,k] <- cECM_vec
        llik[,k] <- llik_vec

      }

    }




  }

 if(inherits(x, "data.frame") ){
   if(!is.null(newdata)){

     ## The model is trained, testing data was included in the call to the cECM function
     return(structure(list(cECM = cECM_out, llik = llik, newdata = newdata.save, rda_fit = rda_fit, x = x, transform = transform), class = "cECM"))

   }else{

     ## The model is trained, but testing data was NOT included in the call to the cECM function
     return(structure(list(rda_fit = rda_fit, x = x, transform = transform), class = "cECM"))

   }
 }else if(inherits(x, "cECM") ){

   if(!is.null(newdata)){

     ## A trained model was provided to the cECM function, as well as testing data.
     return(structure(list(cECM = cECM_out, llik = llik, newdata = newdata.save, rda_fit = rda_fit, x = x$x, transform = x$transform), class = "cECM"))

   }else{

     return(structure(list(rda_fit = rda_fit, x = x$x, transform = x$transform), class = "cECM"))

   }

 }



}

#' Decision Function for the C-ECM Model
#'
#' Returns category decisions for uncategorized events.  When the true category of such events are known performance metrics are returned.
#'
#' @param pval  Class `"cECM"` object obtained from providing training and testing data to the [cECM()] function.
#' @param alphatilde Scalar numeric between 0 and 1, used as the significance level for hypothesis testing.
#' @param vic  Character vector of length one, indicating the Very Important Category (VIC).  Used for judging accuracy, false negatives, and false positives for binary categorization.  Required when `cat_truth` is supplied.
#' @param cat_truth Character vector corresponding to the true group of each row in the `newdata` argument, previously provided to the [cECM()] function.  When supplied, along with `vic`, binary categorization accuracy is returned.
#'
#' @return When `is.null(cat_truth) == TRUE` a vector providing the categorization over all event categories is returned.  When the `cat_truth` and `vic` arguments are supplied a list is returned containing a `data.frame` detailing if each event was categorized accurately in a binary categorization framework, was a false positive, was a false negative, and the estimated event category.  A vector stating overall categorization accuracy, false positive rate, and false negative rate is included with the list.
#'
#' @details
#' When `is.null(cat_truth) == TRUE`, categorization over all event categories is used, using the same framework seen in \insertCite{anderson2007mathematical}{ezECM}.  The return value of `"indeterminant"` happens when there is a failure to reject multiple event categories. `"undefined"` is returned when all event categories are rejected.
#'
#' When the arguments `cat_truth` and `vic` are included, binary categorization is utilized instead of categorization over all training categories.  The definition of accuracy is more ambiguous when categorizing over all training categories and the user is encouraged to develop their own code for such a case.  The goal of binary categorization is to estimate if an uncategorized observation is or is not in the event category stipulated by `vic`.  Uncategorized events which are `"indeterminant"` or `"undefined"` are deemed to not be in the `vic` category.
#'
#'
#' @export
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#'
#' file_path <- system.file("extdata", "good_training.csv", package = "ezECM")
#' training_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
#' newdata <- training_data[1:10,]
#' cat_truth <- newdata$event
#' newdata$event <- NULL
#' training_data <- training_data[-(1:10),]
#'
#' pval <- cECM(training_data, transform = TRUE, newdata = newdata)
#'
#' binary_decision <- cECM_decision(pval = pval, alphatilde = 0.05,
#' vic = "explosion", cat_truth = cat_truth)
#'
#' decision <- cECM_decision(pval = pval, alphatilde = 0.05)
#'
cECM_decision <- function(pval = NULL, alphatilde = NULL, vic = NULL, cat_truth = NULL){

  ### When vic is provided, but cat_truth is not provided, there should be the opportunity to still do binary categorization without checking accuracy metrics

  ## Requirements for the pval argument
  if(!inherits(pval, "cECM")){
    stop("Argument 'pval' must be of class 'cECM' generated by the ezECM::cECM function.")
  }
  if(is.null(pval$cECM)){
    stop("Argument 'pval' does not contain aggregate p-values calculated from uncategorized data.  Use the argument 'newdata' in the ezECM::cECM function to include uncategorized data in the analysis.")
  }

  ## Requirements for the cat_truth and vic arguments

  if(!is.null(cat_truth)){
    if(length(cat_truth) != nrow(pval$cECM)){
      stop("Argument 'cat_truth' must be associated with the data used in previously generating 'pval$cECM'.  It is required that length(cat_truth) == nrow(pval$cECM) ")
    }
    if(is.null(vic)){
      stop("Argument 'vic' is required along with 'cat_truth' in order to calculate accuracy of the categorization in a binary categorization decision framework.")
    }
    if(!is.character(vic)){
      stop("Argument 'vic' must be a character, and must be supplied for calculating accuracy when the 'cat_truth' argument is supplied.")
    }
    if(!(vic %in% names(pval$cECM))){
      stop("Argument 'vic' must be one of the event categories used to train the C-ECM model.")
    }
  }

  if(alphatilde >= 1 | alphatilde <= 0){
    stop("Argument 'alphatilde' must be between the values of 0 and 1, exclusive for both end points.")
  }

  pval <- pval$cECM

  ## Compare the aggregate p-values to alphatilde to determine which categories get rejected

  rejections <- apply(pval, 1, function(X, a){
    names(X)[which(X <= alphatilde)]
  }, a = alphatilde, simplify = FALSE)

  not_rejections <- apply(pval, 1, function(X, a){
    names(X)[which(X > alphatilde)]
  }, a = alphatilde, simplify = FALSE)

  group <- rep(NA, times = length(rejections))

  if(is.null(cat_truth)){

    undefineds <- rep(FALSE, length = length(rejections))
    indeterminants <- rep(FALSE, length = length(rejections))
    ncats <- ncol(pval)

    decision <- rep(NA, times  = length(rejections))

    for(i in 1:length(rejections)){

      if(length(not_rejections[[i]] > 1)){
        ## There is a failure to reject 2 or more event categories
        indeterminants[i] <- TRUE
      }else if(length(rejections[[i]]) == ncats){
        ## All categories are rejected
        undefineds[i] <- TRUE
      }else{
        ## In this case, all but one categories are rejected and this is the decision
        decision[i] <- not_rejections[[i]]
      }

      }

  }else{
    group <- rep(NA, times = length(rejections))
  ## For each uncategorized observation
  for(i in 1:length(rejections)){
    if(vic %in% rejections[[i]]){
      ## very important category is rejected
      group[i] <- "b"
    }else{
      ## very important category is not rejected
      if(length(rejections[[i]]) == (ncol(pval) - 1)){
        ## very important category is the only category not rejected
        group[i] <- "a"
      }else{
        ## very important category, and at least one other category, are not rejected
        group[i] <- "b"
      }
    }
  }
  }



  if(!is.null(cat_truth)){
  true_group <- cat_truth
  true_group[which(cat_truth != vic)] <- "b"
  true_group[which(cat_truth == vic)] <- "a"


    ## return accuracy, false negatives, and false positives, and binary categorization
    out <- data.frame(matrix(NA, ncol = 4, nrow = nrow(pval)))
    names(out) <- c("accuracy", "fn", "fp", "estimated.category")

    out$accuracy <- group == true_group

    out$fn[cat_truth == vic] <- !(out$accuracy[cat_truth == vic])
    out$fp[cat_truth != vic] <- !(out$accuracy[cat_truth != vic])

    out$estimated.category <- group
    out$estimated.category[group == "a"] <- vic
    out$estimated.category[group == "b"] <- paste0("not ", vic)

    metrics <- apply(out[c("accuracy", "fn", "fp")], 2, function(X){
      X <- X[!is.na(X)]
      return(mean(X))
    })

    out <- list(events = out, metrics = metrics)

  }else{

    decision[undefineds] <- "undefined"
    decision[indeterminants] <- "indeterminant"

    out <- decision
  }

  return(out)

}
