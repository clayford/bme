#' Adjusted odds ratio accounting for misclassification
#' 
#' Calculate potential impact of misclassification on an odds ratio estimate for
#' a single 2 x 2 table.
#' 
#' @param data a 2 x 2 matrix with disease=yes in row 1 and exposure=yes in
#'   column 1.
#' @param sensitivity proportion in interval (0,1) of exposed cohort who truly develop the disease
#'   and are diagnosed correctly.
#' @param specificity proportion in interval (0,1) of exposed cohort who truly do not develop the
#'   disease and are diagnosed correctly.
#' @return  A list with two components: 
#'   \describe{ 
#'    \item{Adjusted odds ratio}{Adjusted odds ratio after accounting for misclassification} 
#'    \item{Estimated odds ratio}{Estimated odds ratio based on (possibly) misclassified data} 
#'   }
#' @export
#' @references Newman (2001), pages 72-75, 99.
#' @examples 
#' ## Example 4.2
#' adjusted.odds.ratio(data = breast, sensitivity = 0.9, specificity = 0.99) 
adjusted.odds.ratio <- function(data, sensitivity, specificity){
  if(!is.matrix(data) || nrow(data)!=2 || ncol(data) != 2) stop("data must be a 2 x 2 matrix")
  if(sensitivity > 1 || sensitivity < 0 || specificity > 1 || specificity < 0) stop("sensitivity and specificity must be in (0,1)")
  lhs <- matrix(c(sensitivity, 1-specificity, 1-sensitivity, specificity), nrow = 2, byrow = TRUE)
  rhs <- data[,1]
  low <- solve(lhs, rhs)
  rhs <- data[,2]
  high <- solve(lhs, rhs)
  adj <- (low[1]*high[2])/(low[2]*high[1])
  unadj <- (data[1,1]*data[2,2])/(data[1,2]*data[2,1])
  list("Adjusted odds ratio after accounting for misclassification"=adj,
       "Estimated odds ratio based on (possibly) misclassified data"= unadj)
}

