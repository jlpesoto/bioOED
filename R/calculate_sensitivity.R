
#'
#' Local sensitivities of microbial inactivation
#' 
#' Calculates the local sensitivity function of a microbial inactivation
#' process. These are estimated using finite differences, through the function
#' \code{\link{sensFun}} from the \code{\link{FME}} package.
#' 
#' @param inactivation_model Character defining the inactivation model to use.
#' @param parms Numeric vector with the nominal values of the model parameters.
#' @param n_times Numeric value specifying the nombers of time points where 
#'        the sensitivity functions will be calculated. 100 by default.
#' @param temp_profile Data frame describing the environmental conditions.
#' @param parms_fix Nominal value of the parameters not considered for the
#'        sensitivity.
#' @param varscale The scaling factor for sensitivity variables. \code{NULL}
#' indicates that the variable value is used. 1 by default.
#' @param parscale The scaling factor for parameters. \code{NULL} indicates
#' that the parameter value is used. 1 by default.
#' @param sensvar The output variable for which the sensitivity will be 
#' estimated. \code{"logN"} by default.
#' 
#' @importFrom FME sensFun
#' 
#' @return A data.frame of class \code{sensFun}.
#' 
#' @seealso \code{\link{sensFun}}
#' 
#' @export
#' 
#' @examples
#' parms_fix <- c(temp_ref = 57.5)
#' parms <- c(delta_ref = 3.9,
#'            z = 4.2,
#'            p = 1,
#'            N0 = 1e6
#' )
#' 
#' temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60)
#' )
#' 
#' sensitivity <- sensitivity_inactivation("Mafart", parms,
#'                                temp_profile, parms_fix)
#' 
#' plot(sensitivity)
#' 
sensitivity_inactivation <- function(inactivation_model, parms,
                                     temp_profile, parms_fix,
                                     n_times = 100,
                                     varscale = 1, parscale = 1,
                                     sensvar = "logN", ...) {
    
    times <- seq(0, max(temp_profile$time), length = n_times)
    
    sensitivities <- sensFun(inactivation_sens_handler, parms, 
                             inactivation_model = inactivation_model,
                             times = times, temp_profile = temp_profile,
                             parms_fix = parms_fix, varscale = varscale,
                             parscale = parscale, sensvar = sensvar,
                             ...)
    sensitivities
    
}


#'
#' Handler for the calculation of sensitivities of inactivation models
#' 
#' @importFrom bioinactivation predict_inactivation
#' 
inactivation_sens_handler <- function(model_parms, inactivation_model,
                                      times, temp_profile, parms_fix) { 

    prediction_results <- predict_inactivation(inactivation_model, times,
                                               c(model_parms, parms_fix),
                                               temp_profile)
    
    return(prediction_results$simulation[,1:3])
    
}








