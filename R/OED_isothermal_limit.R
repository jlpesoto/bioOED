
#' Objective function for D-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' 
#' @importFrom dplyr mutate
#' 
detFIM_limit <- function(x, model, pars, limit){
  
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

criterion_Emod <- function(x, model, pars, limit) {
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
criterion_E <- function(x, model, pars, limit) {
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
criterion_Amod <- function(x, model, pars, limit) {
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
criterion_A <- function(x, model, pars, limit) {
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
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' 
#' @inheritParams isothermal_OED
#' 
#' @export
#' 
isothermal_OED_limit <- function(model, pars, limit,
                                 n_points, min_time, max_time, min_temp, max_temp, criterion,
                                 opts = NULL) {
  
  if (min_time <= 0) {
    min_time <- 1e-6
    print("NOTE: min_time has been set to 1e-6 to avoid singularities in Weibullian models")
  }
  if (TRUE) {
    
    if(criterion=="D") {
      
      problem <- list(f = detFIM_limit,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                      ##  x_0 = c(52.53778, 80, 54.13304, 51.1526, 80, 80, 49.63199, 2.526613, 79.78784,
                      ##        10.47441, 52.0859, 52.39074, 60, 52.14292, 52.39072, 52.39071, 52.20695, 60, 52.39645,
                      ##      60 )
      )
      
    }
    if(criterion=="E_mod")
    {
      problem <- list(f = criterion_Emod,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                      ##  x_0 = c(42.32582, 24.96474, 24.87559, 26.01196, 1.089272, 26.03738, 26.13816, 25.8911,
                      ##        21.45262, 24.37624, 53.01773, 52.08012, 52.09364, 52, 59.9312, 52, 52, 52.0151,
                      ##    52.41367, 52.13642   )
      )        
      
    }
    if(criterion=="E")
    {
      problem <- list(f = criterion_E,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                      ## x_0 = c(25.59577, 2.144052, 25.2455, 2.129045, 69.21017, 38.84482, 1.214379, 27.87819, 12.69641,
                      ##       68.77958, 59.85635, 58.77039, 53.43707, 58.78558, 52.88859, 52.50507, 59.99986, 53.22253, 54.92361,
                      ##     53.62613)
      )        
      
    }
    if(criterion=="A_mod")
    {
      problem <- list(f = criterion_Amod,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                      ## x_0 = c(80, 79.83095, 46.50006, 48.1058, 50.54459, 40.30259, 73.10896, 80, 75.6652, 63.59473,
                      ##       52.39071, 52.39528, 53.56418, 53.49075, 53.3838, 53.87354, 52.58552, 52.39071, 52.51119, 52.88706)
      )        
      
    }
    if(criterion=="A")
    {
      problem <- list(f = criterion_A,
                      x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                      x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
                      ##x_0 = c(80, 79.76407, 12.04009, 17.97412, 16.52487, 12.59853, 8.299137, 23.10813, 10.60535, 71.87346, 52.39516, 52.39793, 55.07106,
                      ##      54.16916, 54.37126, 54.9749, 55.88155, 53.62161, 55.34048, 60  )
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
    criterion = "D",
    optim_algorithm = "MEIGO",
    optim_design = my_design,
    limit = limit
  )
  
  class(out) <- c("OEDisothermal", class(out))
  
  out
  
}