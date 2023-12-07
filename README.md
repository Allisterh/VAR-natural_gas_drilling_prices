# README: What Drives Drilling Up and Prices Down?
Replication code and data for 'What Drives Drilling Down and Prices Up? A Structural Vector Autoregressive Model of US Natural Gas Markets'. The code estimates a structural vector autoregressive model of US natural gas markets identified with a mix of sign and zero restrictions using sub-rotations. For the IRFs the modal model (which maximizes the IRF posterior density) as well as 68% highest density bands are plotted. Details can be obtained from the [paper](https://github.com/valwinkler/natural_gas_drilling_prices/blob/main/Winkler_Natural-Gas-Paper.pdf). Feedback for both code and paper is highly appreciated, just [email me](mailto:valentin.winkler@icloud.com).

## R Files

|Category |Filename |Description |
|:------|:------------|:-----------|
|**Main** |`main.R` |Estimate SVAR model, compute statistics/distributions of interest, plot and export|
|**Outsourced** |`bayes_4_4_var.R` |Read in data, estimate SVAR model and compute impulse response functions|
| |`bayes_rotations_4_4.R` |Conduct sub-rotations and check if identifying assumptions are fulfilled given a certain reduced form model|
| |`robustness_functions.R`|Functions for robustness checks |
| |`irfpdf_blockrec.R` |Compute posterior density of structural impulse response function|
|**Auxiliarly** |`fev_decomp.R` |Compute forecast error variance decomposition for SVAR model |
| |`hist_decomp.R` |Compute historical decomposition (conditional on initial values) for SVAR model|
| |`matrix_calc.R` |Useful functions for matrix calculus |
| |`mult_gamma_log.R` |Compute multivariate gamma function (in log-levels)|

## Data Files
|Filename |Description |
|:--------|:-----------|
|`Data_VAR.csv` |Data used in estimating the SVAR model|
|`Data_VAR_Dictionary.csv` |Variable descripiton and sources|
