library(magrittr)

helper_file <- file.path("R", "local-power", "helper-functions-local-power.R")
if (!file.exists(helper_file)) {
	helper_file <- file.path("..", "..", "R", "local-power", "helper-functions-local-power.R")
}

source(helper_file)

# gls_file <- file.path("R", "local-power", "second-stage-gls.R")
# if (!file.exists(gls_file)) {
#   gls_file <- file.path("..", "..", "R", "local-power", "second-stage-gls.R")
# }
# source(gls_file)
