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
  "parallel",          # L
  "parallelly",        # L
  "pbapply",           # L   parallel progress bar
  
  "filelock",          # L   safe concurrent appends from parallel workers
  
  #"abind",              # R   for MFBAVART-SV
  
  "mvtnorm",           # L   for mix$mixBART
  "invgamma"           # L   for mix$mixBART
)

# Pin these versions (required by my code + to prevent stochvol being upgraded by factorstochvol)
pin_versions <- list(
  stochvol       = "3.2.8",        # R\L for sampling stochastic volatilities
  factorstochvol = "1.1.0"         # L   for mix$mixBART , compatible with stochvol
)

# GitHub remotes (package name -> repo)
github_pkgs <- c(
  fatBVARS = "hoanguc3m/fatBVARS"   # L   for BVAR-SV     , "hoanguc3m/fatBVARS" on GitHub
  #mfbvar   = "ankargren/mfbvar"    # R   for MFBAVART-SV , "ankargren/mfbvar" on GitHub
)


# Model switches (set TRUE/FALSE before running)
run_models <- list(
  RW = TRUE,
  RW_SV = TRUE,
  BVAR = FALSE,
  mixBART = FALSE
)

if (!isTRUE(run_models$RW)) {
  stop("RW must be enabled because relative scores are computed versus RW.")
}

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

if (isTRUE(run_models$BVAR)) ensure_github(github_pkgs[["fatBVARS"]], "fatBVARS")

# Re-enforce the critical pin (defensive)
ensure_version("stochvol", pin_versions[["stochvol"]])

# 4) Snapshot once
renv::snapshot(prompt = FALSE)
renv::status()

# 5) Load libraries
pkgs_to_load <- unique(c(
  cran_pkgs,
  names(pin_versions),
  if (isTRUE(run_models$BVAR)) "fatBVARS" else NULL
))
invisible(lapply(pkgs_to_load, library, character.only = TRUE))

set.seed(123) # allow reproducibility

#############################################################################
# MF-BAVART SET-UP (ENVIROMENT and DEPENDENCES)                             #                                                                           
# 1) Path alla repo "vendorizzata" dentro ext/                              #
mf_dir <- here::here("ext", "mf-bavart")                                    #
list.files(mf_dir)
# 2) Crea un environment dedicato per non sporcare il Global Env            #
mf <- new.env(parent = globalenv())                                         #
#
# 3) Carica i file R nell’environment mf                                    #
sys.source(file.path(mf_dir, "aux_func.R"),      envir = mf, chdir = TRUE)  #
sys.source(file.path(mf_dir, "mfbavart_func.R"), envir = mf, chdir = TRUE)  #
#############################################################################

if (isTRUE(run_models$mixBART)) {
  #############################################################################
  # flexBART SET-UP (ENVIROMENT and DEPENDENCES)                              #
  # 1) Path alla repo                                                         #
  mix_dir <- here::here("ext", "flexBART")                                    #
  list.files(mix_dir)
  # 2) Crea un environment dedicato per non sporcare il Global Env            #
  mix <- new.env(parent = globalenv())                                        #
  #
  # 3) Carica i file R nell’environment mix                                   #
  sys.source(file.path(mix_dir, "aux_func.R"), envir = mix, chdir = TRUE)     #
  sys.source(file.path(mix_dir, "flexBART.R"), envir = mix, chdir = TRUE)     #
  #############################################################################
}



###################################################################################
# FUNCTIONS NEEDED (utility helpers for transformations, scoring, and safe I/O):                                                               #
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
# 4) Diebold-Mariano test with Harvey-Leybourne-Newbold small-sample correction.
#    This is used to compare score/loss series of each model against RW.
dm_test_hln <- function(loss_model, loss_rw, h = 1) {
  d <- as.numeric(loss_model) - as.numeric(loss_rw)
  d <- d[is.finite(d)]
  Tn <- length(d)
  if (Tn < 5) {
    return(data.frame(dm_stat = NA_real_, p_value = NA_real_, mean_diff = NA_real_, n_obs = Tn))
  }
  
  d_bar <- mean(d)
  
  gamma0 <- stats::var(d)
  if (!is.finite(gamma0) || gamma0 < 0) gamma0 <- 0
  
  if (h <= 1) {
    var_d_bar <- gamma0 / Tn
  } else {
    max_lag <- min(h - 1, Tn - 1)
    acov_sum <- 0
    if (max_lag > 0) {
      for (lag in seq_len(max_lag)) {
        cov_lag <- stats::cov(d[(lag + 1):Tn], d[1:(Tn - lag)])
        if (is.finite(cov_lag)) acov_sum <- acov_sum + 2 * cov_lag
      }
    }
    long_run_var <- gamma0 + acov_sum
    var_d_bar <- long_run_var / Tn
  }
  
  if (!is.finite(var_d_bar) || var_d_bar <= 0) {
    return(data.frame(dm_stat = NA_real_, p_value = NA_real_, mean_diff = d_bar, n_obs = Tn))
  }
  
  dm_raw <- d_bar / sqrt(var_d_bar)
  hln_factor <- sqrt((Tn + 1 - 2 * h + (h * (h - 1) / Tn)) / Tn)
  dm_stat <- dm_raw * hln_factor
  p_value <- 2 * stats::pt(-abs(dm_stat), df = Tn - 1)
  
  data.frame(dm_stat = dm_stat, p_value = p_value, mean_diff = d_bar, n_obs = Tn)
}
#
# 5) function for Directional Symmetry                                                         #
directional_symmetry <- function(actual, forecast) {                                           #               
  if (length(actual) < 2) {                                                                    #
    return(NA_real_)                                                                           #                                        
  }                                                                                            #                       
  direction_match <- (actual[-1] - actual[-length(actual)]) *                                  #                        
    (forecast[-1] - forecast[-length(forecast)]) > 0                                           #               
  100 / (length(actual) - 1) * sum(direction_match)                                            #              
}                                                                                              #

# 5) Thread-safe CSV appender for parallel workers.
#    Each worker acquires a file lock before appending, so concurrent cycles
#    cannot corrupt files when writing at the same time.
append_rows_to_csv_locked <- function(df, csv_path) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  lock_path <- paste0(csv_path, ".lock")
  lock_obj <- filelock::lock(lock_path, timeout = 60000)
  on.exit(filelock::unlock(lock_obj), add = TRUE)
  
  # Append rows without header (header is created once before parallel loop).
  utils::write.table(
    df,
    file = csv_path,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    append = TRUE
  )
  
  invisible(NULL)
}
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

# Output file paths (written incrementally after each completed cycle)
quantile_scores_path <- here::here("quantile_scores.csv")
qwcrps_scores_path <- here::here("qwcrps_scores.csv")
forecast_history_path <- here::here("forecast_history.csv")

# Recreate output files with headers so workers can safely append rows.
utils::write.table(
  data.frame(tau = numeric(), tau_label = character(), quantile_score = numeric(),
             cycle = integer(), horizon = integer(), model = character()),
  file = quantile_scores_path,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)
utils::write.table(
  data.frame(weighting = character(), qwcrps = numeric(),
             cycle = integer(), horizon = integer(), model = character()),
  file = qwcrps_scores_path,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)
utils::write.table(
  data.frame(cycle = integer(), horizon = integer(), model = character(),
             actual = numeric(), forecast = numeric()),
  file = forecast_history_path,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)


# jx = 222 , 222 is the column representing the 2009M6 data vintage for which at least 228 observations for future contracts data is available.

# 210 row start , and 209 row start when difference variable (we loose one observation otherwise).

run_one_cycle <- function(jx, base_ind) {
  set.seed(100000 + jx)
  cycle_id <- jx - (12 * 18 + 6) + 1
  ind_cycle <- base_ind + (cycle_id - 1)
  
  # Create VAR data in column format 
  rpg <- log(100 * NG_HENRY[210:ind_cycle, jx] / CPI_AC[210:ind_cycle, jx])  # log Real gas price (nominal deflated by US CPI)
  
  for_growth_rpg <- 100 * NG_HENRY[1:ind_cycle, jx] / CPI_AC[1:ind_cycle, jx]  # Real gas price (nominal deflated by US CPI)
  
  # variables in levels 
  capu <- CAPUALF[210:ind_cycle, jx]
  hd <- HDDdev[210:ind_cycle, jx]
  cd <- CDDdev[210:ind_cycle, jx]
  
  # variables in logs
  conslog <- log(CONS_TOT[210:ind_cycle, jx])
  
  # variables in growth rates
  ipg <- 100 * lagn(log(IPALF[1:ind_cycle, jx, drop = FALSE]), 1)
  dryprod <- 100 * lagn(log(PROD_DRY_SA[1:ind_cycle, jx, drop = FALSE]), 1)
  inventories <- 100 * lagn(log(STORE_WORK_SA[1:ind_cycle, jx, drop = FALSE]), 1)
  consg <- 100 * lagn(log(CONS_TOT[1:ind_cycle, jx, drop = FALSE]), 1)
  rigcount <- 100 * lagn(log(RIGS[1:ind_cycle, jx, drop = FALSE]), 1)
  gr_capu <- 100 * lagn(log(CAPUALF[1:ind_cycle, jx, drop = FALSE]), 1)
  gr_rpg <- 100 * lagn(log(for_growth_rpg[1:ind_cycle, jx, drop = FALSE]), 1)
  
  # one observation lost due to differencing
  growth_ip <- ipg[209:nrow(ipg), 1]
  growth_prod <- dryprod[209:nrow(dryprod), 1]
  growth_store <- inventories[209:nrow(inventories), 1]
  growth_cons <- consg[209:nrow(consg), 1]
  growth_rigs <- rigcount[209:nrow(rigcount), 1]
  growth_capu <- gr_capu[209:nrow(gr_capu), 1]
  growth_rpg <- gr_rpg[209:nrow(gr_rpg), 1]
  
  
  # futures contract variable
  futgasrt <- log(fut[210:ind_cycle, ])
  spotgasrt <- log(NG_HENRY[210:ind_cycle, jx])
  inflrt <- log(CPI_AC[(ind_cycle - 120 + 2):ind_cycle, jx]) - log(CPI_AC[(ind_cycle - 120 + 1):(ind_cycle - 1), jx]) # # monthly inflation to calculate average inflation over the past 10 years
  
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
  
  futs <- futgasrt[, 2] - spotgasrt - ((1 + mean(inflrt))^3 - 1) # jj = 2 , we choose to always use future-spot real spread for 3 months contracts always.
  FUTURES_three_months <- futgasrt[, 2]
  
  # Create revised real price of natural gas (most recent vintage)
  x <- 100 * NG_May24[210:(ind_cycle + h), 1] / CPI_May24[210:(ind_cycle + h), 1]  # Real gas price (nominal deflated by US CPI)
  
  # Estimate Models
  data <- list(
    HDD_dev = ts(hd, frequency = 12),                            # monthly
    #DRY_Production = ts(dryprod, frequency = 12),          # monthly
    mgr_DRY_Production = ts(growth_prod, frequency = 12),        # monthly
    #Rig_counts = ts(rigcount, frequency = 12),             # monthly
    mgr_Rig_counts = ts(growth_rigs, frequency = 12),            # monthly
    #CAP_UT = ts(capu, frequency = 12),                     # monthly
    mgr_CAP_UT = ts(growth_capu, frequency = 12),                # monthly
    #Working_Inventories = ts(inventories, frequency = 12), # monthly
    mgr_Working_Inventories = ts(growth_store, frequency = 12),  # monthly
    log_GAS_Price = ts(rpg, frequency = 12)                # monthly
    #mgr_GAS_Price = ts(growth_rpg, frequency = 12)                # monthly
    #CONS = ts(conslog, frequency = 12),                   # monthly
    #mgr_CONS = ts(growth_cons, frequency = 12),                  # monthly
    #INDUSTRIAL_PRODUCTION = ts(ipg, frequency = 12),      # monthly
    #mgr_INDUSTRIAL_PRODUCTION = ts(growth_ip, frequency = 12),   # monthly
    #CDD_dev = ts(cd, frequency = 12),                     # monthly
    #Fut_three = ts(FUTURES_three_months, frequency = 12), # monthly
    #Fut_sp = ts(futs, frequency = 12)                     # monthly
  )
  
  n_draws <- 15000
  t <- length(rpg)
  
  model_draws_by_h <- list()
  
  # RW: random walk without drift and homoskedastic shocks
  if (isTRUE(run_models$RW)) {
    rw_innov <- diff(rpg)
    rw_a0 <- 0.01
    rw_b0 <- 0.01
    rw_a_post <- rw_a0 + length(rw_innov) / 2
    rw_b_post <- rw_b0 + 0.5 * sum(rw_innov^2)
    
    sigma2_draws <- invgamma::rinvgamma(n_draws, shape = rw_a_post, rate = rw_b_post)
    rw_fcst_draws <- matrix(NA_real_, nrow = n_draws, ncol = max(horizons))
    for (i in seq_len(n_draws)) {
      shocks <- rnorm(max(horizons), mean = 0, sd = sqrt(sigma2_draws[i]))
      rw_fcst_draws[i, ] <- rpg[t] + cumsum(shocks)
    }
    model_draws_by_h$RW <- rw_fcst_draws
  }
  
  # RW-SV: random walk with stochastic volatility
  if (isTRUE(run_models$RW_SV)) {
    rwsv_innov <- diff(rpg)
    rwsv_fit <- stochvol::svsample(
      y = rwsv_innov,
      draws = n_draws,
      burnin = 5000,
      quiet = TRUE
    )
    
    latent_h <- rwsv_fit$latent[[1]]
    para <- rwsv_fit$para
    rwsv_fcst_draws <- matrix(NA_real_, nrow = n_draws, ncol = max(horizons))
    
    for (i in seq_len(n_draws)) {
      mu_i <- para[i, "mu"]
      phi_i <- para[i, "phi"]
      sigma_i <- para[i, "sigma"]
      h_prev <- latent_h[i, ncol(latent_h)]
      path <- numeric(max(horizons))
      y_prev <- rpg[t]
      
      for (hh in seq_len(max(horizons))) {
        h_next <- mu_i + phi_i * (h_prev - mu_i) + sigma_i * rnorm(1)
        y_prev <- y_prev + rnorm(1, mean = 0, sd = exp(h_next / 2))
        path[hh] <- y_prev
        h_prev <- h_next
      }
      
      rwsv_fcst_draws[i, ] <- path
    }
    model_draws_by_h$RW_SV <- rwsv_fcst_draws
  }
  
  if (isTRUE(run_models$BVAR)) {
    ##############################################################################
    ### BVAR-SV with Multivariate Skew-t (MST) distribution                   ####
    ##############################################################################
    bvar_data <- do.call(
      cbind,
      lapply(data, function(series) as.numeric(series))
    )
    colnames(bvar_data) <- names(data)
    bvar_data <- bvar_data[complete.cases(bvar_data), , drop = FALSE]
    
    K <- ncol(bvar_data)
    p <- 6
    
    prior <- fatBVARS::get_prior(
      y = bvar_data,
      p = p,
      priorStyle = "Minnesota",
      dist = "MST",
      SV = TRUE,
      lambda1 = 0.2,
      lambda2 = 0.5,
      a_Vprior = 10
    )
    
    inits <- fatBVARS::get_init(prior, samples = n_draws, burnin = 5000)
    
    bvar_fit <- fatBVARS::BVAR.SV(
      y = bvar_data,
      K = K,
      p = p,
      dist = "MST",
      y0 = NULL,
      prior = prior,
      inits = inits
    )
    t_pred <- max(horizons)
    bvar_fcst <- fatBVARS::get_forecast(bvar_fit, t_pred = t_pred, t_current = nrow(bvar_data))
    bvar_fcst_draws <- bvar_fcst$y_pred[horizons, which(colnames(bvar_data) == "log_GAS_Price"), , drop = FALSE]
    
    bvar_matrix <- matrix(NA_real_, nrow = n_draws, ncol = length(horizons))
    for (hh in seq_along(horizons)) {
      bvar_matrix[, hh] <- bvar_fcst_draws[hh, 1, ]
    }
    model_draws_by_h$BVAR <- bvar_matrix
  }
  
  if (isTRUE(run_models$mixBART)) {
    ##############################################################################
    ### mixBART-SV with HS prior                                              ####
    ##############################################################################
    data_mix <- list(
      HDD_dev = ts(hd, frequency = 12),
      #DRY_Production = ts(dryprod, frequency = 12),
      gr_DRY_Production = ts(growth_prod, frequency = 12),
      #Rig_counts = ts(rigcount, frequency = 12),
      gr_Rigs = ts(growth_rigs, frequency = 12),
      #CAP_UT = ts(capu, frequency = 12),
      gr_CAP_UT = ts(growth_capu, frequency = 12),
      #Working_Inventories = ts(inventories, frequency = 12),
      gr_Inventories = ts(growth_store, frequency = 12),
      #Fut_three = ts(FUTURES_three_months, frequency = 12),
      Fut_Spot_sp = ts(futs, frequency = 12),
      #log_GAS_Price = ts(rpg, frequency = 12)
      gr_GAS_Price = ts(growth_rpg, frequency = 12)
    )
    
    flex_data <- do.call(
      cbind,
      lapply(data_mix, function(series) as.numeric(series))
    )
    flex_data <- flex_data[complete.cases(flex_data), , drop = FALSE]
    
    mixbart_fit <- mix$flexBART(
      Yraw = flex_data,
      nburn = 5000,
      nsave = n_draws,
      thinfac = 1,
      prior = "HS",
      prior.sig = c(3, 0.9),
      pr.mean = matrix(0, ncol(flex_data), ncol(flex_data)),
      model = "mixBART",
      sv = "SV",
      fhorz = max(horizons),
      quiet = FALSE
    )
    
    mixbart_fcst_draws <- mixbart_fit$fcst[, horizons, 7, drop = FALSE]
    mixbart_matrix <- matrix(NA_real_, nrow = n_draws, ncol = length(horizons))
    for (hh in seq_along(horizons)) {
      mixbart_matrix[, hh] <- mixbart_fcst_draws[, hh, 1]
    }
    model_draws_by_h$mixBART <- mixbart_matrix
  }
  
  # Evaluate forecasts for requested horizons
  qs_out <- list()
  qw_out <- list()
  fh_out <- list()
  for (h_eval in horizons) {
    h_index <- match(h_eval, horizons)
    actual_level <- x[t + h_eval]
    
    for (model_name in names(model_draws_by_h)) {
      level_draws <- exp(model_draws_by_h[[model_name]][, h_index])
      qs_df <- QS_sample(actual_level, level_draws)
      qs_df$cycle <- cycle_id
      qs_df$horizon <- h_eval
      qs_df$model <- model_name
      qs_out[[length(qs_out) + 1]] <- qs_df
      
      qw_df <- compute_qwcrps_scores(actual_level, level_draws)
      qw_df$cycle <- cycle_id
      qw_df$horizon <- h_eval
      qw_df$model <- model_name
      qw_out[[length(qw_out) + 1]] <- qw_df
      
      point_fcst <- stats::median(level_draws)
      fh_out[[length(fh_out) + 1]] <- data.frame(
        cycle = cycle_id,
        horizon = h_eval,
        model = model_name,
        actual = actual_level,
        forecast = point_fcst
      )
    }
  }
  
  # Combine this cycle outputs and write immediately to disk.
  qs_cycle <- do.call(rbind, qs_out)
  qw_cycle <- do.call(rbind, qw_out)
  fh_cycle <- do.call(rbind, fh_out)
  
  append_rows_to_csv_locked(qs_cycle, quantile_scores_path)
  append_rows_to_csv_locked(qw_cycle, qwcrps_scores_path)
  append_rows_to_csv_locked(fh_cycle, forecast_history_path)
  
  rm(qs_cycle, qw_cycle, fh_cycle, qs_out, qw_out, fh_out, model_draws_by_h)
  gc()
  
  # Return lightweight status only (prevents large in-memory accumulation in master).
  data.frame(cycle = cycle_id, status = "written")
}

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

available_cores <- parallel::detectCores(logical = TRUE)
n_workers <- min(8, available_cores)
cl <- parallel::makeCluster(n_workers, outfile = "")
parallel::clusterSetRNGStream(cl, iseed = 123)
parallel::clusterExport(cl, varlist = "run_models", envir = .GlobalEnv)

parallel::clusterEvalQ(cl, {
  if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
  renv::load()
  
  pkgs_to_load <- unique(c(
    "here", "R.matlab", "dplyr", "pracma", "dbarts", "MASS", "pbapply", "filelock",
    "mvtnorm", "invgamma", "stochvol", "factorstochvol",
    if (isTRUE(run_models$BVAR)) "fatBVARS" else NULL
  ))
  invisible(lapply(pkgs_to_load, library, character.only = TRUE))
  
  if (isTRUE(run_models$mixBART)) {
    mix_dir <- here::here("ext", "flexBART")
    mix <- new.env(parent = globalenv())
    sys.source(file.path(mix_dir, "aux_func.R"), envir = mix, chdir = TRUE)
    sys.source(file.path(mix_dir, "flexBART.R"), envir = mix, chdir = TRUE)
    assign("mix", mix, envir = .GlobalEnv)
  }
  
  NULL
})

base_ind <- ind
jx_vec <- (12 * 18 + 6):(ncol(CPI_AC) - h)

parallel::clusterExport(
  cl,
  varlist = c(
    "NG_HENRY", "CPI_AC", "CAPUALF", "HDDdev", "CDDdev", "CONS_TOT", "IPALF",
    "PROD_DRY_SA", "STORE_WORK_SA", "RIGS", "fut", "NG_May24", "CPI_May24",
    "horizons", "h", "base_ind", "run_one_cycle", "lagn", "QS_sample", "qwCRPS_sample",
    "compute_qwcrps_scores", "quantile_levels", "qwcrps_tau", "qwcrps_weightings", "run_models",
    "append_rows_to_csv_locked", "quantile_scores_path", "qwcrps_scores_path", "forecast_history_path"
  ),
  envir = .GlobalEnv
)

#pbapply::pboptions(type = "timer")
pbapply::pboptions(type = if (interactive()) "timer" else "txt")
res_list <- pbapply::pblapply(jx_vec, run_one_cycle, base_ind = base_ind, cl = cl) # each worker writes its cycle output immediately
#res_list <- parallel::parLapplyLB(cl, jx_vec, run_one_cycle, base_ind = base_ind)
parallel::stopCluster(cl)

# All cycle-level results were already appended to CSV files.
# Read them once for final aggregate evaluation summaries.
quantile_scores <- utils::read.csv(quantile_scores_path)
qwcrps_scores <- utils::read.csv(qwcrps_scores_path)
forecast_history <- utils::read.csv(forecast_history_path)

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

# Diebold-Mariano tests: each model vs RW for each horizon and tau/weighting.
models_vs_rw <- setdiff(unique(quantile_scores$model), "RW")

quantile_dm_results <- list()
for (mdl in models_vs_rw) {
  for (h_eval in sort(unique(quantile_scores$horizon))) {
    for (tau_val in sort(unique(quantile_scores$tau))) {
      rw_loss <- quantile_scores |>
        dplyr::filter(model == "RW", horizon == h_eval, tau == tau_val) |>
        dplyr::arrange(cycle) |>
        dplyr::pull(quantile_score)
      mdl_loss <- quantile_scores |>
        dplyr::filter(model == mdl, horizon == h_eval, tau == tau_val) |>
        dplyr::arrange(cycle) |>
        dplyr::pull(quantile_score)
      
      n_common <- min(length(rw_loss), length(mdl_loss))
      if (n_common == 0) next
      
      dm_out <- dm_test_hln(mdl_loss[seq_len(n_common)], rw_loss[seq_len(n_common)], h = h_eval)
      dm_out$model <- mdl
      dm_out$benchmark <- "RW"
      dm_out$horizon <- h_eval
      dm_out$tau <- tau_val
      quantile_dm_results[[length(quantile_dm_results) + 1]] <- dm_out
    }
  }
}

dm_quantile_scores <- if (length(quantile_dm_results) > 0) do.call(rbind, quantile_dm_results) else data.frame()
utils::write.csv(
  dm_quantile_scores,
  here::here("dm_quantile_scores.csv"),
  row.names = FALSE
)

qwcrps_dm_results <- list()
for (mdl in setdiff(unique(qwcrps_scores$model), "RW")) {
  for (h_eval in sort(unique(qwcrps_scores$horizon))) {
    for (w_name in unique(qwcrps_scores$weighting)) {
      rw_loss <- qwcrps_scores |>
        dplyr::filter(model == "RW", horizon == h_eval, weighting == w_name) |>
        dplyr::arrange(cycle) |>
        dplyr::pull(qwcrps)
      mdl_loss <- qwcrps_scores |>
        dplyr::filter(model == mdl, horizon == h_eval, weighting == w_name) |>
        dplyr::arrange(cycle) |>
        dplyr::pull(qwcrps)
      
      n_common <- min(length(rw_loss), length(mdl_loss))
      if (n_common == 0) next
      
      dm_out <- dm_test_hln(mdl_loss[seq_len(n_common)], rw_loss[seq_len(n_common)], h = h_eval)
      dm_out$model <- mdl
      dm_out$benchmark <- "RW"
      dm_out$horizon <- h_eval
      dm_out$weighting <- w_name
      qwcrps_dm_results[[length(qwcrps_dm_results) + 1]] <- dm_out
    }
  }
}

dm_qwcrps_scores <- if (length(qwcrps_dm_results) > 0) do.call(rbind, qwcrps_dm_results) else data.frame()
utils::write.csv(
  dm_qwcrps_scores,
  here::here("dm_qwcrps_scores.csv"),
  row.names = FALSE
)