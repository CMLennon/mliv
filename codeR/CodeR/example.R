
# The bash script ------------------------------------------------------------------------

# Simulation with simple data (Belloni = F)
Rscript --vanilla run-sim.R F 0

# Simulation with Belloni data (type 1)
Rscript --vanilla run-sim.R T 1

# Simulation with Belloni data (type 2)
Rscript --vanilla run-sim.R T 2

# Simulation with Belloni data (type 3)
Rscript --vanilla run-sim.R T 3

# Logic of the R script ------------------------------------------------------------------

run-sim.R
- Check data type: Belloni T or F
  - F: source("gen-simple.R"); gen_data(type = i[2])
  - T: source("gen-belloni.R"); gen_data(type = i[2])
- If Belloni T: Check beta type (1, 2, or 3)
- Check methods

# The beginning of the R script (run-sim.R) ----------------------------------------------
i = commandArgs(trailingOnly = T)
# Note:
#   - i[1] determines data generation: 'F' for simple and 'T' for Belloni
#   - i[2] gives type of Belloni β distribution
#   - tail(i, -2) gives methods to NOT run
#       - 'rf' for random forest
#       - 'rfcv'
#       - ...

# All methods
methods_all = c("rf", "rfcv", "lasso", "lasso2")
# Methods to run
methods_run = setdiff(methods_all, tolower(tail(i, -2)))

if (i[1] == "F") {
  # Source the function that generates the simple data
  source(here("code", "gen-simple.R"))
  # Define the type of beta distribution (irrelevant for 'simple')
  type = i[2]
} else {
  # Source the function that generates the Belloni data
  source(here("code", "gen-belloni.R"))
  # Define the type of beta distribution
  type = i[2]
}

# Then within iterations:
  # Generate data
  iter_data = gen_data(type)
  # Train model(s)
  ...
