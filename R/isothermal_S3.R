
#' Plotting of IsoSensitivities objects
#' 
#' @importFrom dplyr select 
#' @importFrom tidyr gather %>%
#' @importFrom ggplot2 ggplot geom_raster geom_contour facet_wrap
#' @importFrom ggplot2 scale_fill_gradient2 geom_line
#' 
#' @export
#' 
plot.IsoSensitivities <- function(x, y = NULL, ..., limit = NULL) {
    
    p <- x %>%
        .$sensitivities %>%
        select(times, temperature, ends_with("_scaled")) %>%
        gather(var, sens, -times, -temperature) %>%
        # na.omit() %>%
        ggplot(aes(x = times, y = temperature, z = sens)) +
        geom_raster(aes(fill = sens)) +
        geom_contour(colour = "black") +
        facet_wrap("var") +
        # theme_minimal() +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white")
    
    if (!is.null(limit)) {
        
        p1 <- calculate_limit(x$model, x$pars, limit, range(x$sensitivities$temperature)) %>%
            geom_line(aes(x = tt, y = TT),
                      colour = "green", size = 1, linetype = 2,
                      data = ., inherit.aes = FALSE)
        
        p <- p + p1
        
    }
    
    p
    
}

#' "Detection" limit for each model
#'
#'
calculate_limit <- function(model, pars, limit, temp_range) {
    
    if (model == "Bigelow") {
        
        data.frame(TT = seq(temp_range[1], temp_range[2], length = 100)) %>%
            mutate(tt = limit*pars$D_R*10^(-(.data$TT - pars$temp_ref)/pars$z))
        
    } else if (model == "Peleg") {
        
        data.frame(TT = seq(temp_range[1], temp_range[2], length = 100)) %>%
            mutate(tt = ( limit/log(1 + exp(pars$k_b*(.data$TT - pars$temp_crit)) ) )^(1/pars$n)
                   )
        
    } else if (model == "Mafart") {
        
        data.frame(TT = seq(temp_range[1], temp_range[2], length = 100)) %>%
            mutate(tt = limit^(1/pars$p)*pars$delta_ref*10^(-(.data$TT - pars$temp_ref)/pars$z))
        
    } else {
        stop("Unknown model name")
    }

}

#' Plot of OEDisothermal object
#' 
#' @export
#' 
plot.OEDisothermal <- function(x, y = NULL, ...) {
    
    fake_t <- seq(0, max(x$optim_design$times), length = 50)
    fake_T <- seq(min(x$optim_design$temperature), max(x$optim_design$temperature),
                  length = 50)
    
    p <- expand.grid(fake_t, fake_T) %>%
        rename(times = Var1, temperature = Var2) %>%
        isothermal_sensitivities(x$model, ., x$pars) %>%
        plot(., limit = x$limit)
    
    p + geom_point(aes(times, temperature), inherit.aes = FALSE,
                   data = x$optim_design)
    
}







