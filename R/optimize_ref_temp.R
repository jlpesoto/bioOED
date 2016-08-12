
#'
#' Optimization of the Reference Temperature
#' 
#' Finds the optimum value of the reference temperature which minimizes
#' the correlation between sensitivty functions of the model parameters.
#' 
#' The optimization is made using the \code{\link{optim}} function. The 
#' target for the optimization is the maximization of the determinant
#' of the correlation matrix between parameter sensitivities. The
#' Brent method is used, as it is the recommended one for unidimensional
#' optimization.
#' 
#' @param temp_ref0 Initial value of the reference temperature to use for
#' the optimization.
#' @param lower Lower bound for the reference temperature.
#' @param upper Upper bound for the reference temperature.
#' @param inactivation_model Character identifying the inactivation model
#' to use for the calculation.
#' @param parms Numeric vector with the nominal values of the model parameters.
#' @param n_times Numeric value specifying the nombers of time points where 
#'        the sensitivity functions will be calculated. 100 by default.
#' @param temp_profile Data frame describing the environmental conditions.
#' @param parms_fix Nominal value of the parameters not considered for the
#'        sensitivity.
#'        
#' @return The object returned by \code{\link{optim}}.
#' 
#' @export
#' 
#' @examples
#' 
#' parms <- c(delta_ref = 3.9, z = 4.2, p = 1, N0 = 1e6)
#' temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60))
#' 
#' optimize_refTemp(57, 50, 70, "Mafart", parms, temp_profile, c())
#' 
optimize_refTemp <- function(temp_ref0, lower, upper,
                             inactivation_model, parms,
                             temp_profile, parms_fix,
                             n_times = 100) {
    
    optim(temp_ref0, refTemp_optim_handler,
          lower = lower, upper = upper,
          inactivation_model = inactivation_model,
          parms = parms, temp_profile = temp_profile,
          n_times = n_times,
          parms_fix = parms_fix, method = "Brent",
          control = list(fnscale = -1))
    
    
}

#'
#' Hanlder for the Optimization of Reference Temperature
#' 
refTemp_optim_handler <- function(temp_ref, inactivation_model, parms,
                                  temp_profile, parms_fix, n_times) {
    
    
    parms_fix <- c(temp_ref = temp_ref, parms_fix)
    
    correlations <- calculate_pars_correlation(inactivation_model, parms,
                                               temp_profile, parms_fix,
                                               n_times)
    
    det(correlations)
    
}

#'
#' Correlation Between Model Parameters Sensitivities
#' 
#' @param inactivation_model Character defining the inactivation model to use.
#' @param parms Numeric vector with the nominal values of the model parameters.
#' @param n_times Numeric value specifying the nombers of time points where 
#'        the sensitivity functions will be calculated. 100 by default.
#' @param temp_profile Data frame describing the environmental conditions.
#' @param parms_fix Nominal value of the parameters not considered for the
#'        sensitivity.
#' 
#' @importFrom dplyr select_
#' 
#' @export 
#' 
#' @examples
#' 
#' parms_fix <- c(temp_ref = 57.5)
#' parms <- c(delta_ref = 3.9, z = 4.2, p = 1, N0 = 1e6)
#' temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60))
#' calculate_pars_correlation("Mafart", parms, temp_profile, parms_fix)
#' 
calculate_pars_correlation <- function(inactivation_model, parms,
                                       temp_profile, parms_fix,
                                       n_times = 100) {
    
    sensitivities <- sensitivity_inactivation(inactivation_model, parms,
                                              temp_profile, parms_fix,
                                              n_times)
    
    # sensitivities <- select(sensitivities, -x, -var)
    sensitivities <- select_(sensitivities, quote("-x"), quote("-var"))
    correlations <- cor(sensitivities)
    correlations
}


















