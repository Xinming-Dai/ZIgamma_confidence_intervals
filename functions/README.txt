Read me

Code for "Confidence Interval for Coefficients of Variation of Zero-Inflated Gamma Distribution"

1. "1_bootstrap_method.R", "2_fiducial_method.R", "3_MOVER_method_ver1.R", "4_MOVER_method_ver2.R" have four functions that can construct confidence intervals for the difference of coefficient of variance of two zero-inflated gamma distribution using bootstrap, Fiducial, MOVER1, and MOVER2 methods respectively.

2. simulation_function.R contains `simu` function to calculate coverage probabilities, average lengths, lower errors, and upper errors of the four methods in the mean time. simulation.R `simu` called function in simulation_function.R to simulate.

3. computation_time_simulation.R is for calculating time expense of the four methods.
