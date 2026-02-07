# TESI
MASTER THESIS - APPLIED ECONOMICS and MARKETS (LMAEM) - Nov. 2025

## Description
> This code comes without technical support of any kind. The code is free to use, provided that the paper is cited properly.

The repository contains code files to estimate two different models: 
1) BVAR-SV with Multivariate Skew-t (MST) distribution model (from fatBVARS package);
2) mixBART-SV model (from flexBART(...) function).

Estimation code files for flexBART reside within the folder `flexBART` which resides within the folder `ext`.

### Code files for flexBART
These files create a function flexBART(...) to estimate mixBART-SV model plus other additional models.
The folder `flexBART` contains the estimation functions for different models:
- `flexBART.R` contains the main function flexBART which can estimate the mixBART-SV model.
- `qr.R` contains the main function for bayesian quantile regression model.
The rest of the files are:
- `aux_func.R` collects several auxiliary functions.
- `fcst_script.R` contains an example code for using the function.

### Inputs for flexBART(...)
- `Yraw`            numeric matrix with time observations in rows and variables in columns (T × M).
                  standardized in-place using column means/sds, with na.rm=TRUE allowed, and then used to form the lagged design matrix and response matrix. 
                  This implies:
                  - Rows = time points.
                  - Columns = multiple series.
                  - NAs are allowed (used with na.rm for means/sds).
                  Evidence: standardization and creation of X/Y rely on matrix structure and columnwise means/sds.

- `nburn`           scalar integer (number of burn-in MCMC draws).
                  burn-in length and BART control parameter n.burn.

- `nsave`           scalar integer (number of saved draws).
                  defines total samples ntot, thinning count nthin, and BART control n.samples.

- `thinfac`         (default 1) numeric scalar; thinning multiplier.
                  Used in nthin <- round(thinfac * nsave) to select how many draws are stored.

- `prior`           string flag; the code checks for "Minn" to switch to a Minnesota prior branch.
                  Used as "Minn" triggers Minnesota prior logic for VAR coefficients. Otherwise HS prior machinery is used (default).

- `prior.sig`       numeric vector of length 2 (shape and scale parameters) for the BART residual prior (chisq).
                  Used as chisq(prior.sig[[1]], prior.sig[[2]]) for each equation’s BART residual prior (unless overridden by sv). In sv == "SV" or "heteroBART" it is internally overridden to c(10000^50, 0.5).

- `model`           (default "mixBART") string flag; "mixBART" or "BART" enables BART-driven mean dynamics; other values imply linear VAR-only behavior.
                  Used as branch for whether tree predictions are used when sampling and forecasting.

- `sv`              string flag. Recognized values include "SV", "heteroBART", and "homo" (as seen in script usage).
                  Used as controls stochastic volatility modeling and associated priors; "SV" builds stochastic volatility priors; "heteroBART" uses BART to model volatility. This also changes how prior.sig is set. 

- `fc.approx`       string flag; valid values used in code are "exact" and "approx".
                  Used as switches between exact sampling of the forecast mean (with BART mean + VAR) vs. an approximate form using (A_approx + PHI_draw).

- `restr.var`       NULL or a column name matching colnames(Y); the code uses which(colnames(Y)==restr.var).
                  Used as if non-NULL, forecasts are conditional on specified variable values across quantile grid range.conditional, which adds an extra R dimension to outputs.

- `fhorz`           scalar integer; number of forecast steps.
                  Used as forecast horizon dimension in arrays and loop bounds.

- `quiet`           logical scalar. 
                  Used as if FALSE, shows progress bar and diagnostic plotting during sampling. 

- `pr.mean`         numeric matrix (M × M) used for the prior mean of the linear VAR coefficients.
                  Used as inserted into A_prior for the first M rows/columns, so it must be compatible with the dimension M = ncol(Y).

The function returns a list with two elements:
list("fcst" = fcst_store, "Hfcst" = Hfcst_store)

- `fcst`            numeric array of forecast draws.
                  Dimensions:
                  - If restr.var is NULL: c(nthin, fhorz, M).
                  - If restr.var is not NULL: c(nthin, fhorz, M, R) where R is the number of conditional quantile grid points. 
                  Dimnames:
                  - mcmc1..mcmcN for draw index,
                  - fhorz1..fhorzH for horizon,
                  - variable names from colnames(Y).
                  Scale: stored in the original data scale (the standardization is undone via Ysd and Ymu).
                  Unconditional branch rescales explicitly when storing.
                  Conditional branch stores Y.tp1 * Ysd + Ymu before storing.

- `Hfcst`           numeric array of forecast volatility draws (per series).
                  Dimensions: same as fcst (c(nthin, fhorz, M) or c(nthin, fhorz, M, R)).
                  Dimnames: identical to fcst.
                  Scale: rescaled by Ysd in unconditional branch; conditional branch stores exp(HT) * Ysd.
