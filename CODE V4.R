
# Forecasting natural gas prices in real time

# 0) Make sure renv is available + project library is active
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::load()  # assumes you're inside the renv project

# Optional but recommended if you have a lockfile and want full reproducibility:
if (file.exists("renv.lock")) renv::restore(prompt = FALSE)

# 1) Define requirements
cran_pkgs <- c(
  "here",              # L
  "R.matlab",          # L
  "dplyr",             # L
  "pracma",            # L
  "dbarts",            # R\L sampling trees
  "MASS",              # R\L
  
  "abind",              # R   for MFBAVART-SV
  
  "mvtnorm",           # L   for mix$mixBART
  "invgamma"           # L   for mix$mixBART
)

# Pin these versions (required by my code + to prevent stochvol being upgraded by factorstochvol)
pin_versions <- list(
  stochvol       = "3.2.8",        # 2.0.4", # R\L ver: 2.0.4 for sampling stochastic volatilities
  factorstochvol = "1.1.0"              # 0.9.3"  # L   for mix$mixBART , compatible with stochvol >= 2.0.4
)

# GitHub remotes (package name -> repo)
github_pkgs <- c(
  fatBVARS = "hoanguc3m/fatBVARS", # L   for BVAR-SV     , "hoanguc3m/fatBVARS" on GitHub
  mfbvar   = "ankargren/mfbvar"    # R   for MFBAVART-SV , "ankargren/mfbvar" on GitHub
)


# 2) Helpers
is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

ensure_pkg <- function(pkg) {
  if (!is_installed(pkg)) renv::install(pkg)
}

ensure_version <- function(pkg, ver) {
  if (!is_installed(pkg) || as.character(utils::packageVersion(pkg)) != ver) {
    renv::install(paste0(pkg, "@", ver))
  }
  # hard check
  if (as.character(utils::packageVersion(pkg)) != ver) {
    stop(sprintf("Version pin failed for %s (expected %s).", pkg, ver))
  }
}

ensure_github <- function(repo, pkg_name = NULL) {
  # pkg_name is optional; if omitted we just attempt the install when not present
  if (is.null(pkg_name) || !is_installed(pkg_name)) renv::install(repo)
}

# 3) Install (order matters)
ensure_version("stochvol", pin_versions[["stochvol"]])
ensure_version("factorstochvol", pin_versions[["factorstochvol"]])

for (crnpk in cran_pkgs) ensure_pkg(crnpk)

ensure_github(github_pkgs[["fatBVARS"]], "fatBVARS")
ensure_github(github_pkgs[["mfbvar"]],   "mfbvar")

# Re-enforce the critical pin (defensive)
ensure_version("stochvol", pin_versions[["stochvol"]])

# 4) Snapshot once
renv::snapshot(prompt = FALSE)
renv::status()

# 5) Load libraries
pkgs_to_load <- c(cran_pkgs, names(pin_versions), names(github_pkgs))
invisible(lapply(pkgs_to_load, library, character.only = TRUE))

set.seed(123) # allow reproducibility

#############################################################################
# MF-BAVART SET-UP (ENVIROMENT and DEPENDENCES)                             #                                                                           
# 1) Path alla repo "vendorizzata" dentro ext/                              #
mf_dir <- here::here("ext", "mf-bavart")                                    #
list.files(mf_dir)
# 2) Crea un environment dedicato per non sporcare il Global Env            #                                          #
mf <- new.env(parent = globalenv())                                         #
#
# 3) Carica i file R nell’environment mf                                    #
sys.source(file.path(mf_dir, "aux_func.R"),      envir = mf, chdir = TRUE)  #
sys.source(file.path(mf_dir, "mfbavart_func.R"), envir = mf, chdir = TRUE)  #
#############################################################################

#############################################################################
# flexBART SET-UP (ENVIROMENT and DEPENDENCES)                              #                                                                           
# 1) Path alla repo                                                         #
mix_dir <- here::here("ext", "flexBART")                                    #
list.files(mix_dir)
# 2) Crea un environment dedicato per non sporcare il Global Env            #                                         #
mix <- new.env(parent = globalenv())                                        #
#
# 3) Carica i file R nell’environment mix                                   #
sys.source(file.path(mix_dir, "aux_func.R"), envir = mix, chdir = TRUE)     #
sys.source(file.path(mix_dir, "flexBART.R"), envir = mix, chdir = TRUE)     #
#############################################################################



###################################################################################
# FUNCTIONS NEEDED:                                                               #
# 1) lagn function                                                                #
lagn <- function(data, m) {                                                       #
  # input: data matrix and lag m                                                  #
  data[(m+1):nrow(data), , drop = FALSE] - data[1:(nrow(data)-m), , drop = FALSE] #
}                                                                                 #
#
# 2) function for Quantile Score (QS)                                             #
quantile_levels <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)                       #
QS_sample <-function(true,mcmc,tau=quantile_levels){                              #
  require(pracma)                                                                 #
  #
  tau_len <- length(tau)                                                          #
  Q.tau <- stats::quantile(mcmc,probs=tau)                                        #
  true_rep <- rep(true,tau_len)                                                   #
  QS.vec <- (true_rep-Q.tau)*(tau-((true_rep<=Q.tau)*1))                          #
  #
  data.frame(                                                                     #
    tau       = tau,                                                              #
    tau_label = names(QS.vec),                                                    #
    quantile_score = unname(QS.vec),                                              #
    row.names = NULL                                                              #
  )                                                                               #
}                                                                                 #
#
# 3) function for quantile-weighted CRPS (qwCRPS)                                 #
qwcrps_tau <- seq(0.01, 0.99, by = 0.01)                                          #
qwcrps_weightings <- c("center", "left", "right", "tails")                        #
qwCRPS_sample <-function(true,mcmc,tau=qwcrps_tau,weighting="none"){              #
  require(pracma)                                                                 #
  #
  tau_len <- length(tau)                                                          #
  Q.tau <- stats::quantile(mcmc,probs=tau)                                        #                                       
  true_rep <- rep(true,tau_len)                                                   #                            
  QS.vec <- (true_rep-Q.tau)*(tau-((true_rep<=Q.tau)*1))                          #                                                     
  # 
  weights <- switch(tolower(weighting),                                           #                                      
                    "none" = 1,                                                   #                              
                    "tails" = (2*tau-1)^2,                                        #                                         
                    "right" = tau^2,                                              #                                   
                    "left" = (1-tau)^2,                                           #                                      
                    "center" = tau*(1-tau))                                       #                                        
  wghs <- QS.vec*weights                                                          #                       
  return(pracma::trapz(tau,wghs))                                                 ##############                                 
}                                                                                              #                                                                 
compute_qwcrps_scores <- function(true, mcmc, tau=qwcrps_tau, weightings=qwcrps_weightings) {  #
  data.frame(                                                                                  #
    weighting = weightings,                                                                    #           
    qwcrps = vapply(                                                                           #    
      weightings,                                                                              #    
      function(w) qwCRPS_sample(true, mcmc, tau, weighting = w),                               #                                                
      numeric(1)                                                                               #
    )                                                                                          #
  )                                                                                            #
}                                                                                              #
#
# 4) function for Directional Symmetry                                                         #
directional_symmetry <- function(actual, forecast) {                                           #               
  if (length(actual) < 2) {                                                                    #
    return(NA_real_)                                                                           #                                        
  }                                                                                            #                       
  direction_match <- (actual[-1] - actual[-length(actual)]) *                                  #                        
    (forecast[-1] - forecast[-length(forecast)]) > 0                                           #               
  100 / (length(actual) - 1) * sum(direction_match)                                            #              
}                                                                                              #                     
################################################################################################


###################################################################################

# Load real-time datasets
# Rows: T (starting in 1973M1+36)
# Columns: point in real time (1991M1 to 2024M2)

NG_HENRY <- as.matrix(read.table(here::here("DATA - models", "NG_HENRY.txt"), sep="\t", header=FALSE))	         
# nominal gas price, Henry Hub
CPI_AC <- as.matrix(read.table(here::here("DATA - models", "CPI_AC.txt"), sep="\t", header=FALSE))		           
# average change nowcast of US CPI

# Economic predictor variables:
IPALF <- as.matrix(read.table(here::here("DATA - models", "IPALF.txt"), sep="\t", header=FALSE))		             
# industrial production index (ALFRED vintages)
CAPUALF <- as.matrix(read.table(here::here("DATA - models", "CAPUALF.txt"), sep="\t", header=FALSE))		         
# US capacity utilization rate (ALFRED vintages)

PROD_DRY_SA <- as.matrix(read.table(here::here("DATA - models", "PROD_DRY_SA.txt"), sep="\t", header=FALSE))	   
# marketed NG production
STORE_WORK_SA <- as.matrix(read.table(here::here("DATA - models", "STORE_WORK_SA.txt"), sep="\t", header=FALSE)) 
# underground storage of working gas
CONS_TOT <- as.matrix(read.table(here::here("DATA - models", "CONS_TOT.txt"), sep="\t", header=FALSE))	         
# total NG consumption
RIGS <- as.matrix(read.table(here::here("DATA - models", "RIGS.txt"), sep="\t", header=FALSE))		               
# rig count

CDDdev <- as.matrix(read.table(here::here("DATA - models", "CDDdev.txt"), sep="\t", header=FALSE))		           
# cooling degree days in deviation from historical average
HDDdev <- as.matrix(read.table(here::here("DATA - models", "HDDdev.txt"), sep="\t", header=FALSE))		           
# heating degree days in deviation from historical average

# GAS Futures start in April 1990 (1990M4) but next year gas futures start in June 1990 (1990M6)
gas_futures <- as.matrix(read.table(here::here("DATA - models", "gas_futures.txt"), header=FALSE))
fut <- rbind(matrix(NA, nrow = 207, ncol = 9), gas_futures)

# May 2024 vintage of Henry Hub spot price and CPI (final release data)
# rows: T (1973M1-2024M5)
# columns: most recent vintage of 2024M5
# non-NaN rows start on 1997M1 (row 289)
# C:/Users/HP/Downloads/DATA - LNG/HH_CPI_May2024vintage.mat
HH_CPI <- readMat(here::here("DATA - models", "HH_CPI_May2024vintage.mat"),
                  fixNames = FALSE,          # conserva nomi MATLAB (underscore)
                  drop = "singletonLists",
                  verbose = TRUE)

CPI_May24 <- as.matrix(HH_CPI$CPI_May24)
NG_May24  <- as.matrix(HH_CPI$NG_May24)


# Basic parameters
#ind <- 288+1   # indicates length of initial real-time sample (up to 1997.1)
ind <- 437+1    # indicates length of initial real-time sample (up to 2009.6)
h <- 12         # forecast horizon
horizons <- c(1, 3, 6, 9, 12) # horizons to evaluate
#p <- 6         # fixed lag order

#eval_length <- (ncol(CPI_AC) - h) - (12 * 18 + 6) + 1
#sfe <- matrix(NA, nrow = eval_length, ncol = 2)

quantile_scores <- data.frame()
qwcrps_scores <- data.frame()
forecast_history <- data.frame()

quantile_scores_path <- here::here("quantile_scores.csv")
qwcrps_scores_path <- here::here("qwcrps_scores.csv")
forecast_history_path <- here::here("forecast_history.csv")

#jx = 73   
# 12*6+1 = 73 (6 years and 1 month), is the column representing the 1997M1 data vintage for when we do not use future contracts data.
#jx = 222  
# 222 is the column representing the 2009M6 data vintage for which at least 228 observations for future contracts data is available.
# change 37 with 210 , and change 36 with 209 , when using future contracts data

for (jx in (12*18+6):(ncol(CPI_AC)-h)) {  # adjust for evaluation period 2009M6 to 2024M2
  print(jx)
  cycle_id <- jx - (12 * 18 + 6) + 1
  
  # Create VAR data in column format 
  rpg <- log(100 * NG_HENRY[210:ind, jx] / CPI_AC[210:ind, jx])  # log Real gas price (nominal deflated by US CPI)
  
  # variables in levels 
  capu <- CAPUALF[210:ind, jx]
  hd <- HDDdev[210:ind, jx]
  cd <- CDDdev[210:ind, jx]
  
  # variables in logs
  conslog <- log(CONS_TOT[210:ind, jx])
  
  # variables in growth rates
  ipg <- 100 * lagn(log(IPALF[1:ind, jx, drop = FALSE]), 1)
  dryprod <- 100 * lagn(log(PROD_DRY_SA[1:ind, jx, drop = FALSE]), 1)
  inventories <- 100 * lagn(log(STORE_WORK_SA[1:ind, jx, drop = FALSE]), 1)
  consg <- 100 * lagn(log(CONS_TOT[1:ind, jx, drop = FALSE]), 1)
  rigcount <- 100 * lagn(log(RIGS[1:ind, jx, drop = FALSE]), 1)
  
  
  # one observation lost due to differencing
  ip <- ipg[209:nrow(ipg), 1]
  prod <- dryprod[209:nrow(dryprod), 1]
  store <- inventories[209:nrow(inventories), 1]
  cons <- consg[209:nrow(consg), 1]
  rigs <- rigcount[209:nrow(rigcount), 1]
  
  
  # futures contract variable
  futgasrt <- log(fut[210:ind, ])
  spotgasrt <- log(NG_HENRY[210:ind, jx])
  inflrt <- log(CPI_AC[(ind - 120 + 2):ind, jx]) - log(CPI_AC[(ind - 120 + 1):(ind - 1), jx]) # average inflation over the past 10 years
  #inflrt <- log(CPI_AC[210:ind, jx]) - log(CPI_AC[209:(ind - 1), jx]) # monthly inflation
  
  jj <- switch(as.character(h),
               "1" = 1,
               "3" = 2,
               "6" = 3,
               "9" = 4,
               "12" = 5,
               "15" = 6,
               "18" = 7,
               "21" = 8,
               "24" = 9)
  
  futs <- futgasrt[, 2] - spotgasrt - ((1 + mean(inflrt))^3 - 1)
  #infl_ma <- stats::filter(inflrt, rep(1 / 120, 120), sides = 1) # 10-year rolling average
  #infl_exp <- (1 + infl_ma)^h - 1
  #infl_exp <- (1 + infl_ma)^3 - 1
  #futs <- ts(futgasrt[, jj] - spotgasrt - as.numeric(infl_exp), frequency = 12)
  #futs <- futgasrt[, jj] - spotgasrt - as.numeric(infl_exp)
  #futs <- futgasrt[, 2] - spotgasrt - as.numeric(infl_exp)
  
  # Create revised real price of natural gas (most recent vintage)
  x <- 100 * NG_May24[210:(ind + h), 1] / CPI_May24[210:(ind + h), 1]  # Real gas price (nominal deflated by US CPI)
  
  # Estimate Models
  data <- list(
    log_GAS_Price = ts(rpg, frequency = 12),             # monthly
    mgr_DRY_Production = ts(prod, frequency = 12),       # monthly
    mgr_Working_Inventories = ts(store, frequency = 12), # monthly
    CONS = ts(cons, frequency = 12),                     # monthly
    mgr_Rig_counts = ts(rigs, frequency = 12),           # monthly
    #mgr_INDUSTRIAL_PRODUCTION = ts(ip, frequency = 12),  # monthly
    #CAP_UT = ts(capu, frequency = 12),                   # monthly
    HDD_dev = ts(hd, frequency = 12)                    # monthly
    #CDD_dev = ts(cd, frequency = 12),                    # monthly
    #Fut_sp = ts(futs, frequency = 12)                    # monthly
  )
  
  
  ##############################################################################
  ### BVAR-SV with Multivariate Skew-t (MST) distribution                   ####
  ##############################################################################
  # Estimate BVAR from fatBVARS package BVAR.SV(...)
  bvar_data <- do.call(
    cbind,
    lapply(data, function(series) as.numeric(series))
  )
  colnames(bvar_data) <- names(data)
  bvar_data <- bvar_data[complete.cases(bvar_data), , drop = FALSE]
  
  K <- ncol(bvar_data)
  p <- 6
  h <- h
  
  # Minnesota prior con shrinkage 0.2 e 0.5, MST + SV
  prior <- fatBVARS::get_prior(
    y = bvar_data,
    p = p,
    priorStyle = "Minnesota",
    dist = "MST",
    SV = TRUE,
    lambda1 = 0.2,  # overall shrinkage
    lambda2 = 0.5,  # cross-variable shrinkage
    a_Vprior = 10   # a ~ N(0, 10I)
  )
  
  # Inizializzazione
  inits <- fatBVARS::get_init(prior, samples = 15000, burnin = 5000)
  
  # Stima MST-SV
  bvar_fit <- fatBVARS::BVAR.SV(
    y = bvar_data,
    K = K,
    p = p,
    dist = "MST",
    y0 = NULL,
    prior = prior,
    inits = inits
  )
  t_pred <- max(h, 2)
  bvar_fcst <- get_forecast(bvar_fit, t_pred = t_pred, t_current = nrow(bvar_data))
  bvar_fcst_draws <- bvar_fcst$y_pred[horizons, which(colnames(bvar_data) == "log_GAS_Price"), , drop = FALSE]
  
  
  ##############################################################################
  ### mixBART-SV with HS prior                                              ####
  ##############################################################################
  # Estimate Models
  data <- list(
    log_GAS_Price = ts(rpg, frequency = 12),             # monthly
    mgr_DRY_Production = ts(prod, frequency = 12),       # monthly
    mgr_Working_Inventories = ts(store, frequency = 12), # monthly
    CONS = ts(cons, frequency = 12),                     # monthly
    mgr_Rig_counts = ts(rigs, frequency = 12),           # monthly
    #mgr_INDUSTRIAL_PRODUCTION = ts(ip, frequency = 12), # monthly
    #CAP_UT = ts(capu, frequency = 12),                  # monthly
    HDD_dev = ts(hd, frequency = 12),                    # monthly
    #CDD_dev = ts(cd, frequency = 12),                   # monthly
    Fut_sp = ts(futs, frequency = 12)                    # monthly
  )
  
  # Estimate mixBART from flexBART(...)
  flex_data <- do.call(
    cbind,
    lapply(data, function(series) as.numeric(series))
  )
  flex_data <- flex_data[complete.cases(flex_data), , drop = FALSE] #Y has T= and M=
  
  # fcst.unconditional 
  mixbart_fit <- mix$flexBART(
    Yraw = flex_data,
    nburn = 5000,
    nsave = 15000,
    thinfac = 1,
    prior = "HS",
    prior.sig = c(3, 0.9),
    pr.mean = matrix(0, ncol(flex_data), ncol(flex_data)),
    model = "mixBART",
    sv = "SV",
    #fc.approx="exact",
    fhorz = h,
    quiet = FALSE
  )
  ## This is an array where the first dimension is the number of saved draws, the second is the forecast horizons and the third refer to the different endogenous series
  mixbart_fcst_draws <- mixbart_fit$fcst[, horizons, 1, drop = FALSE]  # draws x horizons x series
  
  
  # Evaluate forecasts for requested horizons
  t <- length(rpg)
  for (h_eval in horizons) {
    h_index <- match(h_eval, horizons)
    actual_level <- x[t + h_eval]
    
    bvar_draws_level <- exp(bvar_fcst_draws[h_index, 1, ])
    mixbart_draws_level <- exp(mixbart_fcst_draws[, h_index, 1])
    rw_draws <- rep(exp(rpg[length(rpg)]), length(bvar_draws_level))
    
    model_draws <- list(
      BVAR = bvar_draws_level,
      mixBART = mixbart_draws_level,
      RW = rw_draws
    )
    
    for (model_name in names(model_draws)) {
      draws <- model_draws[[model_name]]
      qs_df <- QS_sample(actual_level, draws)
      qs_df$cycle <- cycle_id
      qs_df$horizon <- h_eval
      qs_df$model <- model_name
      quantile_scores <- dplyr::bind_rows(quantile_scores, qs_df)
      utils::write.csv(quantile_scores, quantile_scores_path, row.names = FALSE)
      
      qw_df <- compute_qwcrps_scores(actual_level, draws)
      qw_df$cycle <- cycle_id
      qw_df$horizon <- h_eval
      qw_df$model <- model_name
      qwcrps_scores <- dplyr::bind_rows(qwcrps_scores, qw_df)
      utils::write.csv(qwcrps_scores, qwcrps_scores_path, row.names = FALSE)
      
      point_fcst <- stats::median(draws)
      forecast_history <- dplyr::bind_rows(
        forecast_history,
        data.frame(
          cycle = cycle_id,
          horizon = h_eval,
          model = model_name,
          actual = actual_level,
          forecast = point_fcst
        )
      )
      utils::write.csv(forecast_history, forecast_history_path, row.names = FALSE)
    }
  }
  
  # Update index for recursive estimation
  ind <- ind + 1
}

# Evaluate real-time recursive forecast accuracy

# Compute Average Quantile Scores
avg_quantile_scores <- quantile_scores |>
  dplyr::group_by(model, horizon, tau) |>
  dplyr::summarise(mean_score = mean(quantile_score, na.rm = TRUE), .groups = "drop")
rw_quantile_scores <- avg_quantile_scores |>
  dplyr::filter(model == "RW") |>
  dplyr::rename(rw_score = mean_score) |>
  dplyr::select(horizon, tau, rw_score)
avg_quantile_scores_rel <- avg_quantile_scores |>
  dplyr::left_join(rw_quantile_scores, by = c("horizon", "tau")) |>
  dplyr::mutate(relative_to_rw = mean_score / rw_score)
utils::write.csv(
  avg_quantile_scores,
  here::here("avg_quantile_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  avg_quantile_scores_rel,
  here::here("avg_quantile_scores_rel.csv"),
  row.names = FALSE
)

# Compute quantile-weighted Continuous Ranked Probability Scores (qw-CRPS)
avg_qwcrps_scores <- qwcrps_scores |>
  dplyr::group_by(model, horizon, weighting) |>
  dplyr::summarise(mean_score = mean(qwcrps, na.rm = TRUE), .groups = "drop")
rw_qwcrps_scores <- avg_qwcrps_scores |>
  dplyr::filter(model == "RW") |>
  dplyr::rename(rw_score = mean_score) |>
  dplyr::select(horizon, weighting, rw_score)
avg_qwcrps_scores_rel <- avg_qwcrps_scores |>
  dplyr::left_join(rw_qwcrps_scores, by = c("horizon", "weighting")) |>
  dplyr::mutate(relative_to_rw = mean_score / rw_score)
utils::write.csv(
  avg_qwcrps_scores,
  here::here("avg_qwcrps_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  avg_qwcrps_scores_rel,
  here::here("avg_qwcrps_scores_rel.csv"),
  row.names = FALSE
)

# Compute Directional Symmetry (DS) statistic
ds_scores <- forecast_history |>
  dplyr::group_by(model, horizon) |>
  dplyr::summarise(
    ds = directional_symmetry(actual, forecast),
    .groups = "drop"
  )
rw_ds <- ds_scores |>
  dplyr::filter(model == "RW") |>
  dplyr::select(horizon, rw_ds = ds)
ds_scores_rel <- ds_scores |>
  dplyr::left_join(rw_ds, by = "horizon") |>
  dplyr::mutate(relative_to_rw = ds / rw_ds)
utils::write.csv(
  ds_scores,
  here::here("ds_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  ds_scores_rel,
  here::here("ds_scores_rel.csv"),
  row.names = FALSE
)
