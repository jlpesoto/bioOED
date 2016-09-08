
#'
#' Correlation Plot of Parameter Sensitivities
#' 
#' Makes a correlation plot of the sensitivities between model
#' parameters.
#' 
#' @param x Instance of parCorrelation
#' @param y Ignored
#' @param ... Ignored
#' 
#' @importFrom corrplot corrplot
#' 
#' @export
#' 
plot.parCorrelation <- function(x, y = NULL, ...) {
    
    corrplot(x, method="color",  
             type="upper", order="hclust", 
             addCoef.col = "black", # Add coefficient of correlation
             tl.col="black", tl.srt=45, #Text label color and rotation
             # hide correlation coefficient on the principal diagonal
             diag=FALSE )
    
    
    
}









