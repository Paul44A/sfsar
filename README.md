# sfsar
Code and Data for the paper "Spatially Functional Simultaneous Autoregressive Models"

## Folder Structure

- `fos.R`: function-on-scalar baseline estimator
- `fsar_scalar.R`: SFSAR estimator with scalar spatial dependence
- `fsar_func.R`: SFSAR estimator with functional spatial dependence
- `bootstrap_inference.R`: bootstrap inference for Algorithms1--3
- `y_gen.R`: generate spatially dependent response matrix
- `application/`: empirical analysis and data
- `simulation/`: simulation study 

## Empirical application:
- `application/app_inference.R`: empirical estimation and bootstrap inference
- `application/data_1974to2023_3h.RData`: hourly climate data
- `application/sites.RData`: site information
- `application/weightmat_inv.RData`: spatial matrix

## Simulation study:

`simulation/simulation.R`: Estimation accuracy of FoS, SFSAR with scalar rho, SFSAR with rho(t); LM test, F-type test, and pointwise CI for SFSAR with rho(t)

## Required R Packages

- `fda`
- `parallel`
- `Matrix`
- `dplyr`
- `ggplot2`
- `readr`
- `scales`
- `stringr`
- `tidyr`
