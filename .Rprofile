source("renv/activate.R")

if (getwd() == "/home/dvdb/Storage/MesdagProject") {
  source("simulations/rstudio_stan_cmdstanr.R")
  local({

    assignFunctionInPackage <- function(fun, name, package) {
      ns <- getNamespace(package)
      unlockBinding(name, ns)
      assign(name, fun, ns)
      lockBinding(name, ns)
    }

    assignFunctionInPackage(rstudio_stanc_cmdstanr, "rstudio_stanc", "rstan")

  })

}
