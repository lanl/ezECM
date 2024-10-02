## © 2024. Triad National Security, LLC. All rights reserved.
##
## This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
##
## (End of Notice)
##
#' Import p-values
#'
#' Imports and organizes observed p-values located in a `*.csv` file for training or categorization using an existing model.
#'
#' The purpose of this function is to give the user a way to import p-value data in which the data will be organized for use with the [cECM()] and [BayesECM()] functions.  Warnings are printed when potential issues may arise with the supplied file, and the function attempts to detect and correct simple formatting issues.
#'
#' Ideally, the user supplies a `*.csv` file which contains a header row labeling the columns.  Each column contains the p-values of a particular discriminant, and each row must correspond to a single event.  If training data is to be imported, the column containing known event categories is labeled `"event"`.  If new data is imported to be used with an existing model fit, the order of the columns in the `new.data` file must be the same as the `*.csv` file containing training data.
#'
#' @param file Character string providing name of `*.csv` file to be imported.  File name within current working directory or absolute path are acceptable.  Argument passed directly to [utils::read.csv()], see [utils::read.csv()] for details.
#' @param header Logical indicating if first row of supplied `*.csv` file contains column names.  See [utils::read.csv()] for details.
#' @param sep Character string indicating the field separator character for the supplied `*.csv` file.
#' @param training Logical indicating if the supplied `*.csv` file is to be used as training data.  Only serves to suppress warnings.
#'
#' @return A [base::data.frame()] of p-values with each row corresponding to a single event and each column corresponding to a particular discriminant.  If data labels are correctly provided in the supplied `*.csv` file, an additional column labeled `event` will hold these values.
#'
#' @export
#'
#' @importFrom utils read.csv
#'
#' @examples
#'
#' file_path <- system.file("extdata", "good_training.csv", package = "ezECM")
#' training_data <- import_pvals(file = file_path, header = TRUE, sep = ",", training = TRUE)
#'
import_pvals <- function(file = NULL, header = TRUE, sep = ",", training = TRUE){


  dat <- utils::read.csv(file = file, header = header, sep = sep)

  dat[dat == "na"] <- NA

  dat <- as.data.frame(apply(dat, 2, function(X){
    X.numeric <- suppressWarnings(as.numeric(X))
    if(all(is.na(X.numeric))){
      return(X)
    }else{
      return(X.numeric)
    }
  }, simplify = FALSE))

  if(header == TRUE){
    names.split <- strsplit(names(dat), split = "")
    X.detect <- sapply(names.split, function(X){X[1]})
    if(sum(X.detect == "X") == (length(X.detect) - 1)){
      warning("It is likely the first row in the supplied file does not include column names, and it is possible the first row of data was not imported.  Consider labeling columns or setting the argument: header = FALSE")
    }
    if(any(names(dat) == "Event")){
      names(dat)[which(names(dat) == "Event")] <- "event"
    }
    if(length(which(names(dat) == "event")) > 1){
      stop("Multiple columns within .csv file labeled as 'event'.  Only one column must contain event categories.")
    }
    if(length(which(names(dat) == "event")) == 0){
      if(training == TRUE){
        warning("No columns of .csv file labeled as containing event categories.  Imported data cannot be used as training data.  If not the case, ensure the column containing event categories is labeled 'event' (case sensitive).")
      }
    }

  }else{
    event.col <- which(sapply(dat, class) == "character")
    if(length(event.col) >= 2){
      stop("Error in supplied data.  Ensure only one column of the csv file contains character entries, missing data is indicated via NA (case sensitive), and consider supplying a header row and setting the argument: header = TRUE.")
    }else if(length(event.col) == 0){
      if(training == TRUE){
        warning("No column specifying event category detected.  Data cannot be used as training data.  If not the case, consider including a header row in the .csv file and setting the argument: header = TRUE")
      }
    }else{
      names(dat)[event.col] <- "event"
      warning(paste0("Column ", event.col, " assumed to contain event category."))
    }
  }

  if(any(names(dat) == "event")){
    if(sum(is.na(dat$event)) >= 1){
      if(training == TRUE){
        warning("'event' column is missing data.  Supplied .csv data cannot be used to fit an ECM model.")
      }
    }
    dat <- dat[c(names(dat)[-which(names(dat) == "event")], "event")]
  }

return(dat)


}

#' Saving and loading fitted ECM models
#'
#' @description
#'
#' `export_ecm()` saves an ECM model fit to training data for later use. The object can be the output of the [cECM()] or [BayesECM()] functions.  Object saved in user specified location.
#'
#' `import_ecm()` loads a previously exported ECM model into the global environment.
#'
#' @param x Name of ECM model to be exported. `class(x)` should be either `"cECM"` or `"BayesECM"`.
#' @param file Character string containing absolute or relative path of `*.rda` file to be exported or imported.
#' @param verb Logical indicating if a message indicating success should be printed to the console.
#'
#' @return Saves file or imports file into global environment.
#' @export
#'
#' @examples
#'
#' x <- pval_gen(sims = 20, pwave.arrival = list(optim.starts = 5))
#' pval_cat <- cECM(x = x, transform = TRUE)
#'
#' export_ecm(x = pval_cat, file = "example-ecm.rda", verb = FALSE)
#'
#' reload_pval_cat <- import_ecm(file = "example-ecm.rda")
#'
#' ## The below code shouldn't be used,
#' ## it deletes the file created during the example
#'
#' unlink("example-ecm.rda")
#'
#'
#'
export_ecm <- function(x = NULL, file = stop("'file' must be specified"), verb = TRUE){


  save(x, file = file)
  if(verb){
  cat("File: '", file, "' successfully saved.")
  }

}
#' @rdname export_ecm
#' @export
import_ecm <- function(file = NULL, verb = TRUE){

  if(file.exists(file)){
    temp_load <- load(file)
    if(verb){
    message("File successfully loaded")
    }
  }else{
    if(verb){
    message("Cannot find file.  Check name and ensure the correct working directory is being used.")
    }
  }

  return(temp_load)

}
#'
