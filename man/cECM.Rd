% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p-agg.R
\name{cECM}
\alias{cECM}
\title{Multiple Discriminant Analysis}
\usage{
cECM(x, newdata = NULL, rda_params = NULL, transform = TRUE)
}
\arguments{
\item{x}{Either a \code{data.frame} of training data or the \code{list} output from previous use of the \code{cECM} function.  If a \code{data.frame} is supplied, each row contains one point of training data.  Additionally, one column, with the name of \code{event}, must contain labels for the event of each row.  The remaining columns contain the observed \code{ncol(x) - 1} discriminant related p-value data.}

\item{newdata}{a \code{data.frame} containing new data for event categorization.  Must use a subset of the discriminants used in the training data.  If for a particular event a certain discriminant is unavailable, specifying \code{NA} in place of the p-value will allow for calculation of the aggregate p-value.}

\item{rda_params}{a \code{list} of arguments passed to the \code{\link[klaR:rda]{klaR::rda()}} function.  If arguments \code{rda_params$x} or \code{rda_params$grouping} are supplied, they are ignored.}

\item{transform}{Logical indicating if the supplied p-values should be transformed by the function \eqn{2/\pi \times \mathrm{asin}\sqrt{X}}.  Ignored if a \code{list} is supplied as the argument \code{x}.}
}
\value{
A list.  Any returned objects contain a list element indicating the value of \code{transform} supplied to the \code{cECM} function call, as well as a \code{\link[klaR:rda]{klaR::rda()}} object related to relevant training data.  In addition if \code{newdata} argument is supplied, the returned list contains a \code{data.frame} specifying aggregate p-values for each new event (rows) for related event category (columns).
}
\description{
Fits a regularized discriminant analysis model to labeled training data and generates an aggregate p-value for categorizing newly obtained data.
}
\details{
Details on regularized discriminant analysis (RDA) can be found in \insertCite{friedman1989regularized;textual}{ezECM}.  Details on related implementation found in \insertCite{anderson2007mathematical;textual}{ezECM}.
}
\examples{

x <- pval_gen(sims = 20, pwave.arrival = list(optim.starts = 5))
s <- sample(1:20, size = 2)

newdata <- x[s,]
newdata <- newdata[,-which(names(newdata) == "event")]

x <- x[-s,]

pval_cat <- cECM(x = x, transform = TRUE)

pval_cat <- cECM(x = pval_cat, newdata = newdata)


}
\references{
\insertAllCited{}
}
