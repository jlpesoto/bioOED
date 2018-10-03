
#' Objective function for D-optimal OED
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @import tidyselect
#' @importFrom dplyr mutate 
#' @return Numeric value of the objective function for criterium D, which is a determinant of the FIM.
#' 
#' @examples  
#' pars <- list(temp_crit = 55,
#'         n = 1.5,
#'         k_b = 0.1)
#' detFIM(x = c(10,15, 20, 25), "Peleg", pars)
#'
detFIM <- function(x, model, pars){
    
    half <- length(x)/2
    
    time_points <- x[1:half]
    temp_points <- x[(half+1):length(x)]
    
    design <- data.frame(times = time_points, temperature = temp_points)
    -det(calculate_isothermal_FIM(model, design, pars))
}

#' Optimal Experiment Design of isothermal inactivation
#' 
#' OED of microbial inactivation experiments.
#' 
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the nominal model parameters.
#' @param n_points numerical stating the number of data points.
#' @param min_time numerical stating the lower limit for the time points.
#' @param max_time numerical stating the upper limit for the time points.
#' @param min_temp numerical stating the lower limit for the temperature.
#' @param max_temp numerical stating the upper limit for the temperature.
#' @param opts options for the MEIGO algorithm. By default, a maximum of 2000
#' function evaluations with local finish with the DHC algorithm 
#' (see help from MEIGO).
#' 
#' @return A MEIGO object 
#' 
#' 
#' @export
#' 
#' @examples 
#' pars <- list(temp_crit = 55,
#' n = 1.5,
#' k_b = 0.1)
#' OED <- isothermal_OED("Peleg", pars,
#'                             n_points=10, min_time=0, max_time=100, min_temp=52, max_temp=60,
#'                             opts = NULL)
#' OED$optim$xbest
#'
isothermal_OED <- function(model, pars,
                           n_points, min_time, max_time, min_temp, max_temp,
                           opts = NULL) {

    if (min_time <= 0) {
        min_time <- 1e-6
        print("NOTE: min_time has been set to 1e-6 to avoid singularities in Weibullian models")
    }
    if (TRUE) {
        
        problem <- list(f = detFIM,
                        x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                        x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                        )
    }
    
    if (is.null(opts)) {
        
        opts <- list(maxeval=2000,local_finish="DHC")
    }
    
    result <- MEIGO(problem, opts, algorithm="ESS", model = model, pars = pars)
    
    ## Build the design matrix
    
    half <- length(result$xbest)/2
    
    time_points <- result$xbest[1:half]
    temp_points <- result$xbest[(half+1):length(result$xbest)]
    
    my_design <- data.frame(times = time_points, temperature = temp_points)
    
    ## Return
    
    out <- list(
        optim = result,
        model = model,
        pars = pars,
        criteria = "D",
        optim_algorithm = "MEIGO",
        optim_design = my_design,
        limit = NULL
    )
    
    class(out) <- c("OEDisothermal", class(out))
    
    out
}























