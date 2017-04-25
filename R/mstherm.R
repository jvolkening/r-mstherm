#' @title Model and analyze MS/MS-based protein melting data.
#'
#' @description \code{mstherm} is a package for modeling and analysis of
#' MS/MS-based thermal proteome profiling (TPP) experiments.
#'
#' @name mstherm
#' @author Jeremy Volkening \email{jdv@@base2bio.com}
#' @docType package
#'
#' @importFrom stats coefficients density dnorm loess mad median nls.control
#'  p.adjust pnorm quantile var lm
#' @importFrom graphics par plot points rect mtext lines curve legend arrows
#'  abline
#' @importFrom utils read.delim txtProgressBar setTxtProgressBar
#'  getTxtProgressBar
#' @importFrom parallel detectCores
#' @importFrom foreach foreach "%dopar%"
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom RColorBrewer brewer.pal
#' @importFrom nls2 nls2
#' @importFrom plotrix addtable2plot
NULL
