#' Simulation of Landcover change over time, based on sedimentary pollen record time series
#'
#' ....
#'
#' @name LandCoverChange-package
#' @docType package
#' @author Peter Ewald and Simeon Lisovski.
NULL

## Data import

##' Read data
##'
##' The \code{readWhatever} function is and example function.
##' @title Read Whatever
##' @param file the light file to import.
##' @param skip number of initial lines to skip
##' @return Returns a dataframe with columns
##' \item{\code{data.frame}}{a valid data frame if all goes well}
##' @importFrom utils read.csv
##' @export
readWhatever <- function(file, skip = 0) {
  d <- read.csv(file,header=FALSE,skip=skip)
  d
}