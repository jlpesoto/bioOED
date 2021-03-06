% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OED_penalty.R
\name{inactivation_OED_penalty}
\alias{inactivation_OED_penalty}
\title{Optimum Experimental Design of Microbial Inactivation with Penalty}
\usage{
inactivation_OED_penalty(inactivation_model, parms, temp_profile,
  parms_fix, n_points, time_min, criteria = "D", n_times = 100,
  sensvar = "logN", optim_algorithm = "global", opts_global = NULL,
  ...)
}
\arguments{
\item{inactivation_model}{Character string defining the inacivation model.}

\item{parms}{Named numeric vector defining the model parameters. They must
be named according to the needs of \code{\link{predict_inactivation}}.}

\item{temp_profile}{Data frame defining the temperature profile. It must 
contain a column named \code{time} and a column named
\code{temperature}.}

\item{parms_fix}{Named numeric vector defining the model parameters to be 
omitted during the calculation of the local sensitivities.}

\item{n_points}{Number of measurements which will be taken during the
experiment.}

\item{time_min}{Numeric value indicating the minimum space between
measurements.}

\item{criteria}{Character defining the criteria for the OED. Either 
\code{D} (default) or \code{E-mod}.}

\item{n_times}{Integer defining th enumber of discrete time points used for
the interpolation of the local sensitivities.}

\item{sensvar}{Character defining the variable to use for the OED. Either
\code{logN} (default) or \code{N}.}

\item{optim_algorithm}{Character defining the type of algorithm to use for
the optimization. Either \code{global} (default) or \code{local}.}

\item{opts_global}{List defining the options for the global optimization
algorithm (see \code{\link{MEIGO}}). By default, global solver with
a maximum of 50000 function evaluations and printout on every step.}

\item{...}{Additional arguments passed to penalty_function.}
}
\value{
A list of class \code{OEDinactivation} with the following items:
     \itemize{
         \item optim: Objetc returned by the optimization function.
         \item model: Inactivation model used for the calculations.
         \item parms: Nominal model parameters.
         \item parms_fix: Model parameters not considered for the
               sensitivity calculation.
         \item criteria: Criteria used for the OED.
         \item sensvar: Variable used for the OED.
         \item optim_algorithm: Type of optimization algorithm.
         \item optim_times: Optimum measurement times calculated.
         \item penalty: Logical indicating whether penalty function was
               used.
         \item temp_profile: Temperature profile of the experiment.
         }
}
\description{
Performs an optimum experimental design for the settings selected including
a function which penalties points too close. The
OED is based on the FIM, estimated using the local sensitivity functions
provided by \code{\link{sensitivity_inactivation}}.
}
\examples{

## Definition of input variables

parms_fix <- c(temp_ref = 57.5)
parms <- c(delta_ref = 3.9,
           z = 4.2,
           p = 1,
           N0 = 1e6
)

temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60))

n_points <- 5
time_min <- 10

## OED with local optimization

set.seed(0123182)

local_OED <- inactivation_OED_penalty("Mafart", parms, temp_profile, parms_fix,
                              n_points, criteria = "E-mod", sensvar = "logN",
                              optim_algorithm = "local", time_min = time_min)

print(local_OED$optim_times)
plot(local_OED)

## OED with global optimization

opts_global <- list(maxeval=500,  local_solver=0,
                    local_finish="DHC", local_iterprint=1)

global_OED <- inactivation_OED_penalty("Mafart", parms, temp_profile, parms_fix,
                               n_points, criteria = "E-mod", opts_global = opts_global,
                               time_min = time_min)

print(global_OED$optim_times)
plot(global_OED)


}
