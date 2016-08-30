
#'
#' Optimum Experimental Design of Microbial Inactivation
#' 
#' @importFrom MEIGOR MEIGO
#' 
#' @export
#'
make_OED <- function(inactivation_model, parms, temp_profile, parms_fix,
                     n_points, criteria = "D",
                     n_times = 100, sensvar = "logN",
                     optim_algorithm = "local") {
    
    ## Calculate sensitivities
    
    sensitivities <- sensitivity_inactivation(inactivation_model, parms,
                                              temp_profile, parms_fix,
                                              n_times = n_times, sensvar = sensvar)
    
    ## Prepare for optimization
    
    if (grepl(criteria, "D")) {
        
        tgt_function <- objective_D
        
    } else if (grepl(criteria, "E-mod")) {
        
        tgt_function <- objective_Emod
        
    } else {
        
        stop(paste("Unknown criteria:", criteria))
        
    }
    
    max_time <- max(temp_profile$temp)
    problem <- list(f=tgt_function,
                    x_L=rep(1e-3, n_points),  # 1e-3 to avoid singularity
                    x_U=rep(max_time, n_points)
                    )
    
    opts <- list(maxeval=50000,  local_solver=0,
                 local_finish="DHC", local_iterprint=1)
    
    ## Make optimization
    
    # times <- seq(1, max_time, length = n_points)
    times <- runif(n_points, 1, max_time)
    
    if (grepl(optim_algorithm, "local")) {
        
        results <- optim(times, tgt_function,
                         lower = problem$x_L, upper = problem$x_U,
                         sensitivities = sensitivities,
                         method = "L-BFGS-B")
        
    } else {
        results <- MEIGO(problem, opts, algorithm="ESS",
                         sensitivities = sensitivities)
        
    }

    
    ## Return results
    
    results

    
}

#'
#' Objective Function for the D Criterium
#' 
objective_D <- function(times, sensitivities) {
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- criterium_D(FIM)
    
    out
    
}

#'
#' Objective Function for the modified-E Criterium
#' 
objective_Emod <- function(times, sensitivities) {
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- criterium_modE(FIM)
    
    out
    
}








