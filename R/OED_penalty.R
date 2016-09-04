
#'
#' Penalty Function for OED
#' 
penalty_function <- function(time_points, time_min, a = 1e15, b = 2e15) {
    
    sorted_times <- sort(time_points)
    differences <- diff(sorted_times)
    min_diff <- min(differences)
    
    if (time_min > min_diff) {
        
        results <- ((a-b)/3)*min_diff + b
        
    } 
    else {
        
        results <- 0
        
    }
    
    results
    
}

#'
#' Objective Function for the D Criterium with Penalty
#' 
objective_D_penalty <- function(times, sensitivities, time_min) {
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- criterium_D(FIM) + penalty_function(times, time_min)
    
    out
    
}

#'
#' Objective Function for the modified-E Criterium with Penalty
#' 
objective_Emod_penalty <- function(times, sensitivities, time_min) {
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- criterium_modE(FIM) + penalty_function(times, time_min)
    
    out
    
}

#'
#' Optimum Experimental Design of Microbial Inactivation with Penalty
#' 
#' @export
#' 
inactivation_OED_penalty <- function(inactivation_model, parms, temp_profile, parms_fix,
                                     n_points, time_min, criteria = "D",
                                     n_times = 100, sensvar = "logN",
                                     optim_algorithm = "global",
                                     opts_global = NULL) {
    
    ## Calculate sensitivities
    
    sensitivities <- sensitivity_inactivation(inactivation_model, parms,
                                              temp_profile, parms_fix,
                                              n_times = n_times, sensvar = sensvar)
    
    ## Prepare for optimization
    
    if (grepl(criteria, "D")) {
        
        tgt_function <- objective_D_penalty
        
    } else if (grepl(criteria, "E-mod")) {
        
        tgt_function <- objective_Emod_penalty
        
    } else {
        
        stop(paste("Unknown criteria:", criteria))
        
    }
    
    max_time <- max(temp_profile$temp)
    problem <- list(f=tgt_function,
                    x_L=rep(1e-3, n_points),  # 1e-3 to avoid singularity
                    x_U=rep(max_time, n_points)
    )
    
    if (is.null(opts_global)) {
        opts_global <- list(maxeval=50000,  local_solver=0,
                            local_finish="DHC", local_iterprint=1)
    }
    
    ## Make optimization
    
    if (grepl(optim_algorithm, "local")) {
        
        times <- runif(n_points, 0, max_time)
        
        results <- optim(times, tgt_function,
                         lower = problem$x_L, upper = problem$x_U,
                         sensitivities = sensitivities,
                         method = "L-BFGS-B", time_min = time_min)
        
        best_points <- results$par
        
    } else if (grepl(optim_algorithm, "global")) {
        results <- MEIGO(problem, opts_global, algorithm="ESS",
                         sensitivities = sensitivities, time_min = time_min)
        
        best_points <- results$xbest
        
    } else {
        stop(paste("Unknown optimization algorithm:", optim_algorithm))
    }
    
    ## Return results
    
    out <- list(optim = results)
    
    class(out) <- c("OEDinactivation", class(out))
    
    out$model <- inactivation_model
    out$parms <- parms
    out$parms_fix <- parms_fix
    out$criteria <- criteria
    out$sensvar <- sensvar
    out$optim_algorithm <- optim_algorithm
    out$optim_times <- sort(best_points)
    out$penalty <- TRUE
    out$time_min <- time_min
    
    out
}






















