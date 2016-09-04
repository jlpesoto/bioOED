
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


