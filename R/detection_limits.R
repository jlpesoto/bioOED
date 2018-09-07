

#' Detection limit of the Bigelow model
#' 
detection_bigelow <- function(pars, temperature, limit) {
    
    limit * pars$D_R * 10^(-(temperature - pars$temp_ref)/pars$z) 
}

#' Detection limit of the Peleg model
#' 
detection_peleg <- function(pars, temperature, limit) {
    
    ( limit/log(1 + exp(pars$k_b*(temperature - pars$temp_crit)) ) )^(1/pars$n)
    
}

#' Detection limit of the Mafart model
#'
#'
detection_mafart <- function(pars, temperature, limit) {
    
    limit^(1/pars$p) * pars$delta_ref * 10^(-(temperature - pars$temp_ref)/pars$z) 
}

#' Calculate detection limit
#' 
get_detection <- function(model, pars, temperature, limit) {
    
    switch(model,
           Bigelow = detection_bigelow(pars, temperature, limit),
           Mafart = detection_mafart(pars, temperature, limit),
           Peleg = detection_peleg(pars, temperature, limit)
           )
}

# calculate_limit("Peleg", pars, 6, c(52, 60)) %>%
#     mutate(aa = get_detection("Peleg", pars, TT, 6))







