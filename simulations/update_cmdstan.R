# This file is only used to update cmdstan

cpp_options <- list(
  "CXXFLAGS+= -oFast -march=native -mtune=native"
)
cmdstanr::install_cmdstan(cpp_options = cpp_options, cores = 24)
cmdstanr::cmdstan_make_local(cpp_options = cpp_options)
cmdstanr::rebuild_cmdstan(cores = 24)
