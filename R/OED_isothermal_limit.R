
#'
#' Objective function for D-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @import tidyselect
#' @importFrom dplyr mutate 
#' @return Numeric value of the objective function for criterium D, which is a determinant of the FIM.
#' 
#' @examples  
#' pars <- list(temp_crit = 55,
#'         n = 1.5,
#'         k_b = 0.1)
#' criterium_D(x = c(10,15, 20, 25), "Peleg", pars, limit=7)
#'
criterium_D_iso <- function(x, model, pars, limit){
  
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select(times, temperature)
  
  -det(calculate_isothermal_FIM(model, design, pars))
}

#'
#' Objective function for E modified-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium E modified, which is a determinant of the FIM.
#' @examples  
#' pars <- list(temp_crit = 55,
#'         n = 1.5,
#'         k_b = 0.1)
#' criterium_Emod(x = c(10,15, 20, 25), "Peleg", pars, limit=7)
#' 
criterium_Emod_iso <- function(x, model, pars, limit) {
  tol_eigen <- 1e-100
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select(times, temperature)
  eigenvalues <- eigen(calculate_isothermal_FIM(model, design, pars))
  if(abs(min(eigenvalues$values))-tol_eigen<0){
    return(1e6)
  }
  else if(abs(min(eigenvalues$values))-tol_eigen>0) {
    return(abs(max(eigenvalues$values)/min(eigenvalues$values)))
  }
}

#'
#' Objective function for E-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium E, which is a determinant of the FIM.
#'
#' @examples  
#' pars <- list(temp_crit = 55,
#'         n = 1.5,
#'         k_b = 0.1)
#' criterium_E(x = c(10,15, 20, 25), "Peleg", pars, limit=7)
#' 
criterium_E_iso <- function(x, model, pars, limit) {
  tol_det <- 1e-5
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select(times, temperature)
  if(abs(det(calculate_isothermal_FIM(model, design, pars)))-tol_det<0) {
    return(1e100)
  }
  else {
    eigenvalues <- eigen(solve(calculate_isothermal_FIM(model, design, pars)))
    return(max(eigenvalues$values))
    
  }
}

#'
#' Objective function for A modified-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium A modified, which is a determinant of the FIM.
#' 
#' @examples  
#' pars <- list(temp_crit = 55,
#'         n = 1.5,
#'         k_b = 0.1)
#' criterium_Amod(x = c(10,15, 20, 25), "Peleg", pars, limit=7)
#'
criterium_Amod_iso <- function(x, model, pars, limit) {
  half <- length(x)/2
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select(times, temperature)
  return(-sum(diag(calculate_isothermal_FIM(model, design, pars))))
  
}

#'
#' Objective function for A-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium A, which is a determinant of the FIM.
#' 
#' @examples  
#' pars <- list(temp_crit = 55,
#'         n = 1.5,
#'         k_b = 0.1)
#' criterium_A(x = c(10,15, 20, 25), "Peleg", pars, limit=7)
#'
criterium_A_iso <- function(x, model, pars, limit) {
  tol_det <- 1e-5
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select(times, temperature)
  if(abs(det(calculate_isothermal_FIM(model, design, pars)))-tol_det<0) {
    return(1e100)
  }
  else {
    return(sum(diag(solve(calculate_isothermal_FIM(model, design, pars)))))
    
  }
}
#'
#' OED of isothermal microbial inactivation with detection limit
#' 
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the nominal model parameters.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @param n_points numerical stating the number of data points.
#' @param min_time numerical stating the lower limit for the time points.
#' @param max_time numerical stating the upper limit for the time points.
#' @param min_temp numerical stating the lower limit for the temperature.
#' @param max_temp numerical stating the upper limit for the temperature.
#' @param criterium character string defining the criterium to use.
#' @param opts options for the MEIGO algorithm. By default, a maximum of 2000
#' function evaluations with local finish with the DHC algorithm 
#' (see help from MEIGO).
#' 
#' @return A MEIGO object 
#' 
#' @import tidyverse
#' 
#' @export
#' 
#' @examples 
#' pars <- list(temp_crit = 55,
#' n = 1.5,
#' k_b = 0.1)
#' OED <- isothermal_OED_limit("Peleg", pars, limit=7,
#'                             n_points=10, min_time=0, max_time=100, min_temp=52, max_temp=60, criterium="D",
#'                             opts = NULL)
#' OED$optim$xbest
#'

isothermal_OED_limit <- function(model, pars, limit,
                                 n_points, min_time, max_time, min_temp, max_temp, criterium,
                                 opts = NULL) {
  
  if (min_time <= 0) {
    min_time <- 1e-6
    print("NOTE: min_time has been set to 1e-6 to avoid singularities in Weibullian models")
  }
  if (TRUE) {
    
    if(criterium=="D") {
      
      problem <- list(f = criterium_D_iso,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                     
      )
      
    }
    if(criterium=="E_mod")
    {
      problem <- list(f = criterium_Emod_iso,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))

      )        
      
    }
    if(criterium=="E")
    {
      problem <- list(f = criterium_E_iso,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))

      )        
      
    }
    if(criterium=="A_mod")
    {
      problem <- list(f = criterium_Amod_iso,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                      
      )        
      
    }
    if(criterium=="A")
    {
      problem <- list(f = criterium_A_iso,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))

      )        
      
    }
  }
  
  if (is.null(opts)) {
    
    opts <- list(maxeval=2000,local_finish="DHC")
  }
  
  result <- MEIGO(problem, opts, algorithm="ESS",
                  model = model, pars = pars, limit = limit)
  
  ## Map the results back
  
  half <- length(result$xbest)/2
  
  time_points <- result$xbest[1:half]
  temp_points <- result$xbest[(half+1):length(result$xbest)]
  
  my_design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select(times, temperature)
  
  ## Return
  
  out <- list(
    optim = result,
    model = model,
    pars = pars,
    criterium = "D",
    optim_algorithm = "MEIGO",
    optim_design = my_design,
    limit = limit
  )
  
  class(out) <- c("OEDisothermal", class(out))
  
  out
  
}