

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
                                 n_points, min_time, max_time, min_temp, max_temp,
                                 opts = NULL) {
    
    if (min_time <= 0) {
        min_time <- 1e-6
        print("NOTE: min_time has been set to 1e-6 to avoid singularities in Weibullian models")
    }
    if (TRUE) {
        
        problem <- list(f = detFIM_limit,
                        x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                        x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
        )
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
        criteria = "D",
        optim_algorithm = "MEIGO",
        optim_design = my_design,
        limit = limit
    )
    
    class(out) <- c("OEDisothermal", class(out))
    
    out
    
}
