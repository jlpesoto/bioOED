
#'
#' Calculation of Fisher Information Matrix
#' 
#' @param sensitivities data.frame of class \code{sensFun} as returned by
#' \code{\link{sensitivity_inactivation}}.
#' 
#' @return Matrix with the estimation of the Fisher Information Matrix.
#' 
#' @importFrom dplyr select_
#' @importFrom dplyr %>%
#' 
#' @export
#' 
calculate_FIM <- function(sensitivities) {
    
    sensitivities <- select_(sensitivities, quote("-x"), quote("-var")) %>%
        as.matrix()
    
    FIM <- t(sensitivities) %*% sensitivities
    FIM
}



















