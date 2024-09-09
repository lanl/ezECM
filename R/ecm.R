## © 2024. Triad National Security, LLC. All rights reserved.
##
## This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
##
## (End of Notice)
##
#' Plot the data and output of [cECM()] categorization
#'
#' @param x an object of which is of class `"cECM"`, retrieved from the output of the [cECM()] function.  The `"cECM"` object may or may not contain unlabeled data.
#' @param discriminants character or integer vector of length two.  If a character vector is provided, the character strings must match a subset of the column names for the training data, ie. `all(discriminants %in% names(x$x))` must be true.  If an integer vector is given, the elements of the vector select the column indices of the data to be plotted.
#' @param thenull character string or `NULL`.  When unlabeled data is found within an `"cECM"` object, the name of one of the event categories can be provided as this argument.  When `"thenull"` is provided, unlabled data where the category hypothesis is rejected is colored red.
#' @param alphatilde numeric value specifying hypothesis testing significance level.  Used in conjunction with `thenull`, aggregate p-values less than `alphatilde` are rejected and colored accordingly.
#' @param ... arguments passed to [base::plot()]
#'
#' @details
#'
#' The plot generated from plot.ecm() is first dependent on if the provided `x$newdata` contains a data frame of unlabled data.
#'
#' If unlabled data is not part of the `"cECM"` object, the labled data is simply plotted with the 0.68 and 0.95 confidence levels obtained from the distribution fits returned from [cECM()].
#'
#' If unlabled data ***is*** part of the `"cECM"` object, the unlabled data is plotted in addition to the distribution plots.  Each unlabled data point appears on the plot as an integer, which indexes the corresponding row of `x$newdata`.
#'
#' @return Plot illustrating results of [cECM()]
#' @export
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
#' plot(x = pval_cat, thenull = "explosion")
#'
#' @importFrom grDevices hcl.colors
#' @importFrom graphics points lines par legend text
#' @importFrom ellipse ellipse
#' @importFrom stats cor cov2cor
#'
#' @method plot cECM
#'
plot.cECM <- function(x, discriminants = c(1,2), thenull = NULL, alphatilde = 0.05,...){


    alphatilde <- 1-alphatilde

    ## Ensuring par settings are not changed upon exit
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))

    if(length(discriminants) != 2){
      stop("Please provide a vector of length 2 for the discriminant argument.")
    }

    if(is.character(discriminants)){
      if(!(all(discriminants %in% names(x$x)))){
        stop("Discriminant names not found within training data set labels.")
      }

      if(!is.null(x$newdata)){
        if(!all(discriminants %in% names(x$newdata))){
          stop("Discriminant names not found within newdata.  Make sure organization of newdata matches training data and consider refitting the model.")
        }
        newdata_sub <- x$newdata[discriminants]
      }
      training_sub <- x$x[discriminants]
    }else if(is.numeric(discriminants)){
      if(any(names(x$x)[discriminants] == "event")){
        stop("Invalid option:  'event' column chosen for plotting.  Change 'dicriminants' argument to only subset discriminants.")
      }

      if(!is.null(x$newdata)){
        newdata_sub <- x$newdata[, discriminants]
      }

      training_sub <- x$x[,discriminants]

      discriminants <- names(training_sub)

    }else{
      stop("Please provide a numeric or character vector for the discriminant argument.")
    }


      plot.args <- list(...)
      plot.args <- as.list(plot.args)

      if(is.null(plot.args$main)){
        if(x$transform == TRUE){
          plot.args$main <- "p-values transformed"
        }else{
          plot.args$main <- "p-values"
        }
      }

      if(is.null(plot.args$xlab)){
        plot.args$xlab <- discriminants[1]
      }

      if(is.null(plot.args$ylab)){
        plot.args$ylab <- discriminants[2]
      }

      if(is.null(plot.args$ylim)){
        plot.args$ylim <- c(0,1)
      }

      if(is.null(plot.args$xlim)){
        plot.args$xlim <- c(0,1)
      }

      x.cats <- list()

      for(i in unique(x$x$event)){
        x.cats[[i]] <- x$x[which(x$x$event == i),]
      }

      cols.use <- grDevices::hcl.colors(n = length(unique(x$x$event)), palette = "Zissou 1")

      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = FALSE)

      do.call(plot, c(1, type = "n", plot.args))

      for(i in unique(x$x$event)){
        # plot ellipses related to fit
        Sig <- x$rda_fit$covariances[discriminants ,discriminants ,i]
        mu <- x$rda_fit$means[discriminants, i]

        ell68 <- ellipse::ellipse(stats::cov2cor(Sig), scale = sqrt(diag(Sig)), centre = mu,
                                  level = c(0.68))

        ell95 <- ellipse::ellipse(stats::cov2cor(Sig), scale = sqrt(diag(Sig)), centre = mu,
                                  level = c(0.95))

        graphics::lines(ell68[,1], ell68[,2], lty = 3, col = cols.use[which(unique(x$x$event) == i)])
        graphics::lines(ell95[,1], ell95[,2], lty = 1, col = cols.use[which(unique(x$x$event) == i)])

      }

      for(i in unique(x$x$event)){
        # plot points
        x.cat <- x$x[which(x$x$event == i), discriminants]
        graphics::points(x.cat[,1], x.cat[,2], pch = 19, col = cols.use[which(unique(x$x$event) == i)])
      }

      if(!is.null(x$newdata)){
        newdata.filled <- x$newdata[, discriminants]
        missing.data <- data.frame(t(apply(x$newdata[,discriminants],1,function(X){is.na(X)})))

        if(is.null(thenull)){
          ## If null hypothesis is not specified, simply plot the points at their locations
          newdata.filled[is.na(newdata.filled)] <- 0.5

          points(newdata.filled[missing.data[,1], discriminants[1]], newdata.filled[missing.data[,1], discriminants[2]], pch = 6, lwd= 0.5, cex = 2.5)
          points(newdata.filled[missing.data[,2], discriminants[1]], newdata.filled[missing.data[,2], discriminants[2]], pch = 2, lwd= 0.5, cex = 2.5)

          graphics::text(newdata.filled[, discriminants[1]], newdata.filled[, discriminants[2]], labels = (1:nrow(x$newdata)), adj = c(0.5, 0.5))


          graphics::par(xpd = TRUE)
          graphics::legend("topright", inset=c(-0.28,0), legend= c(unique(x$x$event), "New Data", paste0(discriminants[1], "\nMissing"), paste0(discriminants[2], "\nMissing")),
                           pch=c(rep(19, times = length(unique(x$x$event))),49, 6, 2), col = c(cols.use, "black", "black", "black"),
                           pt.cex = c(rep(1, times = length(unique(x$x$event)) + 1), 2.5, 2.5), y.intersp = 1.8)
        }else{
          ## If null hypothesis is specified, plot if rejected.

          if(!(thenull %in% unique(x$x$event))){
            stop("If provided, argument 'thenull' must match event categories specified in training data.")
          }

          for(i in discriminants){
            newdata.filled[is.na(newdata.filled[,i]),i] <- x$rda_fit$means[which(row.names(x$rda_fit$means) == i), which(colnames(x$rda_fit$means) == thenull)]
          }

          failed <- which(x$cECM[, thenull] < 1-alphatilde)

          if(length(failed) != nrow(x$cECM)){
          graphics::text(newdata.filled[-failed, discriminants[1]], newdata.filled[-failed, discriminants[2]], labels = (1:nrow(x$newdata))[-failed], adj = c(0.5, 0.5))
            }
          if(length(failed) != 0){
          graphics::text(newdata.filled[failed, discriminants[1]], newdata.filled[failed, discriminants[2]], labels = (1:nrow(x$newdata))[failed], adj = c(0.5, 0.5), col = "red")
            }
          graphics::par(xpd = TRUE)
          legend.loc <- graphics::legend("topright", inset=c(-0.32,0), legend= c(unique(x$x$event), paste0("H0: ", thenull),   "Rejected", "Failed to\nReject"),
                           pch=c(rep(19, times = length(unique(x$x$event))), NA, 49, 50), col = c(cols.use, "black", "red", "black"),
                           x.intersp = c(rep(1, times = length(unique(x$x$event))), -0.5, 1,1),
                           y.intersp = c(rep(1, times = length(unique(x$x$event)) + 2),2), plot = TRUE)

          lines(x = legend.loc$rect$left + c(0, legend.loc$rec$w * 0.8),
                y = rep(legend.loc$text$y[(length(unique(x$x$event)) + 1):(length(unique(x$x$event)) + 2)] %*% c(0.72, 0.28), times = 2), lty = 3)
        }
      }else{
        graphics::legend("topright", inset=c(-0.3,0), legend= unique(x$x$event),
                         pch=19, col = cols.use, xpd = TRUE)
      }


}
