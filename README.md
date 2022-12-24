# ZIgamma_confidence_intervals

This repository contains simulation functions, simulation results, and plots for the paper—Confidence Intervals for Coefficients of Variance of Zero-inflated Gamma Distribution (Xinming Dai).

## Abstract

In metrology, variability of precipitation is often demonstrated utilizing coefficient of variation. Furthermore, the distinction of rainfall from two areas or times can be demonstrated by the difference of two independent coefficients of variation. Rainfall data is often right-skewed and has massive zeros, and tests show rainfall data is distributed as zero-inflated gamma. Common existing methods of inferencing the confidence interval for the zero-inflated gamma distribution is unsatisfactory because the information regarding zeros is disregarded. In this work, the fiducial generalized confidence interval, the parametric bootstrap, and MOVER are proposed to solve the confidence interval for the difference of two coefficients of variance of zero-inflated gamma distribution, and Monte Carlo simulation gives an insight into gauging these methods. By avoiding complicated math computation, these methods are easy to calculate. Monte Carlo simulation indicates that the fiducial method outperforms others in terms of balancing coverage probability and average length. Bootstrap confidence intervals perform well for maintaining balanced tail error rates. In terms of time efficiency, fiducial and MOVER are conspicuously more efficient than bootstrap method as sample size gets large. All methods are exemplified utilizing data involving 11 years monthly rainfall in Beijing and Zhengzhou distributed as zero-inflated gamma distribution.
