# model_uncertainty
Code for estimating model structure errors in memory-carrying time series models assuming the error models of observational data are known.
Case study for the Bates Catchment

See inside scripts folder for code for running the method in the Adaptive Metropolis scheme.
state_uncertainty_trial_gr4j_sma_likelihood_prior_AM.r contains the prior and likelihood functions of the method.

The main script is:
Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood_master.r 
which runs:
Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood_runner.r
which in turn runs:
Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood.r

Input data is contained in the data folder

The rainfall-runoff model code is in the packages folder.
