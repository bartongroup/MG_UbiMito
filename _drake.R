source("R/packages.R")
source("R/data.R")
source("R/compare.R")
source("drake/plan.R")


plan <- bind_rows(
  get_data,
  get_numbers,
  compare_data
)

cfg <- drake_config(plan)