.PHONY: simulation
	
simulation: results/simulations/intermediate-objects/p_values_tbl.rds
	
results/simulations/intermediate-objects/p_values_tbl.rds: R/simulations/simulations.R results/simulations/intermediate-objects/scenarios_dgm_tbl.rds R/simulations/simulation-A4LEARN-setup.R
	Rscript R/simulations/simulations.R
	
results/simulations/intermediate-objects/scenarios_dgm_tbl.rds: R/simulations/simulation-A4LEARN-setup.R R/simulations/data-generating-mechanism.R
	Rscript R/simulations/data-generating-mechanism.R

	
	