# Define required packages
required_packages <- c(
  "dplyr", "truncnorm", "MASS", "ggplot2", "survival", "survminer",
  "conflicted", "rstan", "coda", "cmdstanr", "posterior", "bayesplot",
  "stringr", "ggrepel", "tidyverse", "ggpubr", "loo", "kableExtra",
  "tableone", "gridExtra"
)

# Install any missing packages
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set conflict preferences
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")


# ---- Tumor Size Simulation Functions ----

# Function to simulate longitudinal tumor size data with random effects and error
simulate_sld <- function(n, mu_M0, mu_kg, mu_ks, omega_M0, omega_kg, omega_ks, beta, Z, 
                         time_points_per_patient, error_sd = 0.07, M0_cap = 220, seed_num = 4) {
  # Generate random effects with truncated normal distribution
  set.seed(seed_num)
  random_effects <- matrix(NA, n, 3)
  for (i in 1:3) {
    random_effects[, i] <- rnorm(n, mean = 0, sd = 1) * c(omega_M0, omega_kg, omega_ks)[i]
  }  
  
  # Initialize lists for storing individual patient data
  sld_data <- list()
  M0_data <- list()
  kg_data <- list()
  ks_data <- list()
  
  
  # Loop through each patient to generate tumor size over time
  for (i in 1:n) {
    M0_i <- min(exp(log(mu_M0) + random_effects[i, 1]), M0_cap)  # Baseline tumor size with cap
    kg_i <- exp(log(mu_kg) + random_effects[i, 2])               # Growth rate
    ks_i <- exp(log(mu_ks) + random_effects[i, 3] + beta * Z[i]) # Shrinkage rate
    
    # Initialize a data frame to store SLD for each time point for this patient
    time_points <- time_points_per_patient[[i]]
    sld <- data.frame(id = i, time = time_points, sld = rep(NA, length(time_points)), Z = Z[i])
    
    set.seed(seed_num * i)
    # Simulate SLD for follow-up time points with proportional error
    for (j in 1:length(time_points)) {
      t <- time_points[j]
      if (t < 0) {
        sld$sld[j] <- M0_i * exp(kg_i * t) # Tumor growth before treatment
      } else {
        sld$sld[j] <- M0_i * (exp(kg_i * t) + exp(-ks_i * t) - 1) # Tumor growth + shrinkage
      }
      
      # Add proportional error to the SLD measurement
      sld$sld[j] <- max(sld$sld[j] * (1 + rnorm(1, mean = 0, sd = error_sd)), 0)
    }
    
    # Append patient data to list
    sld_data[[i]] <- sld
    M0_data[[i]] <- data.frame(M0_i = rep(M0_i, length(sld$sld)))
    ks_data[[i]] <- data.frame(ks_i = rep(ks_i, length(sld$sld)))
    kg_data[[i]] <- data.frame(kg_i = rep(kg_i, length(sld$sld)))
  }
  
  # Combine all patients' data into one data frame
  sld_data_df <- do.call(rbind, sld_data)
  M0_data_df <- do.call(rbind, M0_data)
  kg_data_df <- do.call(rbind, kg_data)
  ks_data_df <- do.call(rbind, ks_data)
  
  data_df <- do.call(cbind, list(sld_data_df, M0_data_df, kg_data_df, ks_data_df))
  return(data_df)
}


# ---- Tumor Size Filtering Function ----

# Function to apply the 20% rule and additional stochastic stopping
apply_sld_rule <- function(sld_data, stop_prob = 0.07, seed_num = 123) {
  set.seed(seed_num)
  filtered_data <- sld_data %>%
    group_by(id) %>%
    arrange(time) %>%
    mutate(
      # Track the minimum SLD after time 0
      min_sld_after_time0 = cummin(ifelse(time > 0, sld, Inf)),
      # Calculate the percentage difference from the smallest SLD after time 0
      sld_diff = ifelse(time > 0, (sld / min_sld_after_time0) - 1, 0),
      # Flag when increase > 20% from minimum (this will be the last observed measurement)
      sld_remove_flag = ifelse(time > 0 & sld_diff > 0.2, 1, 0),
      # Create a flag to indicate the first time the 20% rule is triggered
      sld_remove_flag = cummax(sld_remove_flag),  # Cumulative max to propagate the flag forward
      # Introduce a random stopping flag with a given probability at each time point
      last_obs_flag = cumsum(sld_remove_flag) == 1,
      # Introduce a random stopping flag with a given probability at each time point
      random_stop = ifelse(row_number() > 2, rbinom(n = 1, size = 1, prob = stop_prob), 0),
      # Flag all time points after a random stop occurs
      random_stop_flag = cumsum(random_stop)
    ) %>%
    mutate(random_stop_flag = ifelse(random_stop_flag > 1, 1, random_stop_flag)) |> 
    filter(random_stop_flag != 1) |> 
    # Combine the 20% rule (last_obs_flag) and the random stopping flag
    rowwise() |> 
    mutate(final_remove_flag = sum(sld_remove_flag, random_stop_flag)) %>%
    # Retain only the first time a 20% increase is observed (last_obs_flag == TRUE)
    mutate(keep_last_obs = ifelse(last_obs_flag, 1, 0)) %>%
    rowwise() |> 
    mutate(final_remove_flag_2 = final_remove_flag - keep_last_obs) %>%
    # Filter out rows where the final flag indicates stopping/removal but keep the last observed measurement
    filter(final_remove_flag_2 == 0 ) %>%
    # Remove unnecessary columns
    dplyr::select(
      -sld_diff, -sld_remove_flag, -min_sld_after_time0, -random_stop,
      -random_stop_flag, -final_remove_flag, -last_obs_flag, -keep_last_obs
    )
  
  return(filtered_data)
}

# ---- Tumor Growth Model Derivatives ----

# Calculate time to nadir
time_to_nadir <- function(k_g, k_s, beta, Z) {
  t_nadir <- log(k_s / k_g) / (k_g + k_s)
  return(t_nadir)
}

# Calculate the rate of change in tumor size (slope) over time
sld_slope <- function(k_g, k_s, M0, t, beta, Z) {
  slope <- M0 * (k_g * exp(k_g * t) - k_s * exp(-k_s * t))
  return(slope)
}

# ---- Joint Model Survival Simulation ----

# Weibull survival simulation function
simulate_survival_jm_weibull <- function(n, lambda, kappa, gamma, beta_os, Z, time_points_per_patient,
                                         ks_params, kg_params, M0_values, link_function, beta, seed_num = 4) {
  survival_times <- rep(NA, n)
  event_observed <- rep(1, n) # No censoring, all events are observed
  set.seed(seed_num)
  for (i in 1:n) {
    time_points <- time_points_per_patient[[i]]
    # Compute the hazard depending on the selected link function
    if (link_function == "current_sld") {
      # Calculate SLD dynamically based on the tumor model at each time point
      link_values <- sapply(time_points, function(t) {
        if (t < 0) {
          # Tumor growth before treatment
          return(M0_values[i] * exp(kg_params[i] * t))
        } else {
          # Tumor growth + shrinkage after treatment
          return(M0_values[i] * (exp(kg_params[i] * t) + exp(-(ks_params[i]) * t) - 1))
        }
      })
    } else if (link_function == "sld_slope") {
      # Current SLD slope as link function
      link_values <- sapply(time_points, function(t) sld_slope(kg_params[i], ks_params[i], M0_values[i], t))
    } else if (link_function == "ks_parameter") {
      # Shrinkage rate ks as link function (constant over time)
      link_values <- rep(log(ks_params[i]), length(time_points))
    } else if (link_function == "time_to_nadir") {
      # Time to nadir as link function (constant over time)
      t_nadir <- time_to_nadir(kg_params[i], ks_params[i])
      link_values <- rep(t_nadir, length(time_points))
    } else {
      stop("Unknown link function.")
    }
    # Similar code for other link functions...
    cumulative_hazard <- sapply(time_points[time_points > 0], function(t) {
      H_t <- (t / lambda)^kappa
      H_t <- H_t * exp(gamma * link_values[which(time_points == t)] + beta_os * Z[i])
      return(H_t)
    })
    survival_prob <- exp(-cumulative_hazard)
    event_time <- min(time_points[time_points > 0][survival_prob < runif(n = 1, min = 0, max = 1)], na.rm = TRUE)
    survival_times[i] <- event_time
    if (survival_times[i] > 2) event_observed[i] <- 0
    if (survival_times[i] > 2) survival_times[i] <- 2
  }
  return(list(survival_times = survival_times, event_observed = event_observed))
}

generate_block_randomization <- function(n_patients, block_size = 8) {
  # Generate balanced treatment assignments with block randomization
  n_blocks <- ceiling(n_patients / block_size)
  Z <- rep(c(rep(0, block_size / 2), rep(1, block_size / 2)), n_blocks)
  Z <- sample(Z)[1:n_patients]  # Shuffle to randomize, then subset to n_patients
  
  # Ensure the first patient is assigned to the control group
  Z[1] <- 0
  return(Z)
}



# ---- Summary Statistics Function ----

# Function to calculate summary statistics for tumor size and survival data
generate_summary_stats <- function(sld_data, survival_data) {
  # Calculate the number of patients in each group
  patients_control <- nrow(survival_data %>% filter(Z == 0) %>% distinct(id))
  patients_treatment <- nrow(survival_data %>% filter(Z == 1) %>% distinct(id))
  
  # Calculate the number of SLD observations in each group
  sld_measurements_control <- nrow(sld_data %>% filter(Z == 0))
  sld_measurements_treatment <- nrow(sld_data %>% filter(Z == 1))
  
  # Calculate the min, median, and max SLD in each group
  min_sld_control <- min(sld_data %>% filter(Z == 0) %>% pull(sld))
  max_sld_control <- max(sld_data %>% filter(Z == 0) %>% pull(sld))
  min_sld_treatment <- min(sld_data %>% filter(Z == 1) %>% pull(sld))
  max_sld_treatment <- max(sld_data %>% filter(Z == 1) %>% pull(sld))
  
  # Calculate median and range of SLD observations per patient in each group
  median_sld_obs_control <- sld_data %>%
    filter(Z == 0) %>%
    group_by(id) %>%
    tally() %>%
    summarise(median = median(n), min = min(n), max = max(n))
  
  median_sld_obs_treatment <- sld_data %>%
    filter(Z == 1) %>%
    group_by(id) %>%
    tally() %>%
    summarise(median = median(n), min = min(n), max = max(n))
  
  # Median survival time
  median_survival_control <- summary(survfit(Surv(time = survival_data$survival_times, event = survival_data$status) ~ Z, data = survival_data))$table["Z=0", "median"]
  median_survival_treatment <- summary(survfit(Surv(time = survival_data$survival_times, event = survival_data$status) ~ Z, data = survival_data))$table["Z=1", "median"]
  
  # Count of events and censored patients in each group
  events_control <- nrow(survival_data %>% filter(Z == 0 & status == 1))
  censored_control <- nrow(survival_data %>% filter(Z == 0 & status == 0))
  events_treatment <- nrow(survival_data %>% filter(Z == 1 & status == 1))
  censored_treatment <- nrow(survival_data %>% filter(Z == 1 & status == 0))
  
  # Combine all statistics into a list
  summary_stats <- list(
    patients_control = patients_control,
    patients_treatment = patients_treatment,
    sld_measurements_control = sld_measurements_control,
    sld_measurements_treatment = sld_measurements_treatment,
    min_sld_control = round(min_sld_control),
    max_sld_control = round(max_sld_control),
    min_sld_treatment = round(min_sld_treatment),
    max_sld_treatment = round(max_sld_treatment),
    median_sld_obs_control = median_sld_obs_control,
    median_sld_obs_treatment = median_sld_obs_treatment,
    median_survival_control = median_survival_control,
    median_survival_treatment = median_survival_treatment,
    events_control = events_control,
    censored_control = censored_control,
    events_treatment = events_treatment,
    censored_treatment = censored_treatment
  )
  
  return(summary_stats)
}

# ---- LaTeX Table Generator Function ----

# Function to generate a LaTeX-formatted table based on the summary statistics
generate_latex_table <- function(summary_stats, scenario_name, sample_size) {
  # Extract values from the summary statistics list
  latex_code <- paste0(
    "\\begin{table}[h]\n",
    "    \\centering\n",
    "    \\begin{tabular}{lcc}\n",
    "        \\toprule\n",
    "        \\textbf{Metric} & \\textbf{Control Group} & \\textbf{Treatment Group} \\\\\n",
    "        \\midrule\n",
    "        Number of patients & ", summary_stats$patients_control, " & ", summary_stats$patients_treatment, " \\\\\n",
    "        Number of SLD observations & ", summary_stats$sld_measurements_control, " & ", summary_stats$sld_measurements_treatment, " \\\\\n",
    "        SLD range & ", summary_stats$min_sld_control, "-", summary_stats$max_sld_control, " & ", summary_stats$min_sld_treatment, "-", summary_stats$max_sld_treatment, " \\\\\n",
    "        Median number of SLD observations per patient & ", 
    summary_stats$median_sld_obs_control$median, " (", summary_stats$median_sld_obs_control$min, "-", summary_stats$median_sld_obs_control$max, ") & ",
    summary_stats$median_sld_obs_treatment$median, " (", summary_stats$median_sld_obs_treatment$min, "-", summary_stats$median_sld_obs_treatment$max, ") \\\\\n",
    "        Median survival time (years) & ", summary_stats$median_survival_control, " & ", summary_stats$median_survival_treatment, " \\\\\n",
    "        Number of events & ", summary_stats$events_control, " & ", summary_stats$events_treatment, " \\\\\n",
    "        Number of censored patients & ", summary_stats$censored_control, " & ", summary_stats$censored_treatment, " \\\\\n",
    "        \\bottomrule\n",
    "    \\end{tabular}\n",
    "    \\caption{Summary Statistics for ", scenario_name, ", Sample Size: ", sample_size, ".}\n",
    "\\end{table}\n"
  )
  
  return(latex_code)
}



# Function to calculate individual survival probabilities over time with all ks and Z combinations
calculate_individual_survival <- function(time_points, M0, k_g, k_s, Z, lambda, kappa, gamma, beta_os, link_function, beta) {
  
  # Adjust k_s based on treatment status
  if (Z == 1) {
    ks_control <- k_s * exp(-beta)     # Control group shrinkage rate for treated individuals
    ks_experimental <- k_s             # Experimental group shrinkage rate for treated individuals
  } else {
    ks_control <- k_s                  # Control group shrinkage rate for control individuals
    ks_experimental <- k_s * exp(beta) # Experimental group shrinkage rate for control individuals
  }
  
  # Define a nested function to compute the link value based on k_s values
  get_link_value <- function(t, k_s, k_g, M0, link_function = "current_sld") {
    if (link_function == "current_sld") {
      if (t < 0) {
        return(M0 * exp(k_g * t))  # Before treatment: only growth
      } else {
        return(M0 * (exp(k_g * t) + exp(-k_s * t) - 1))  # After treatment: growth + shrinkage
      }
      
    } else if (link_function == "sld_slope") {
      return(M0 * (k_g * exp(k_g * t) - k_s * exp(-k_s * t)))
      
    } else if (link_function == "ks_parameter") {
      return(log(k_s))
      
    } else if (link_function == "time_to_nadir") {
      t_nadir <- log(k_s / k_g) / (k_g + k_s)
      return(t_nadir)
      
    } else {
      stop("Unknown link function specified.")
    }
  }
  
  # Calculate survival probabilities based on different ks and Z values
  # calculate_survival_prob <- function(current_k_s, current_Z) {
  #   sapply(time_points, function(t) {
  #     # Cumulative hazard up to time t
  #     link_value <- get_link_value(t, k_s, k_g, M0)
  #     H_t <- (t / lambda)^kappa
  #     H_t <- H_t * exp(gamma * link_value + beta_os * current_Z)
  #     # Survival probability as exponent of negative cumulative hazard
  #     exp(-H_t)
  #   })
  # }
  
  # Survival probabilities for all combinations
  survival_prob_control_Z0 <- calculate_survival_prob(ks_control, 0)
  survival_prob_control_Z1 <- calculate_survival_prob(ks_control, 1)
  survival_prob_experimental_Z0 <- calculate_survival_prob(ks_experimental, 0)
  survival_prob_experimental_Z1 <- calculate_survival_prob(ks_experimental, 1)
  
  
  NIE <- survival_prob_experimental_Z1 - survival_prob_control_Z1
  NDE <- survival_prob_control_Z1 - survival_prob_control_Z0
  TE <- survival_prob_experimental_Z1 - survival_prob_control_Z0
  PTE <- ifelse(survival_prob_experimental_Z1 != survival_prob_control_Z0, 
                (survival_prob_experimental_Z1 - survival_prob_control_Z1) / (survival_prob_experimental_Z1 - survival_prob_control_Z0), 
                NA)
  
  
  # Return survival probabilities in a data frame
  return(data.frame(
    time = time_points,
    survival_prob_control_Z0 = survival_prob_control_Z0,
    survival_prob_control_Z1 = survival_prob_control_Z1,
    survival_prob_experimental_Z0 = survival_prob_experimental_Z0,
    survival_prob_experimental_Z1 = survival_prob_experimental_Z1,
    NIE = NIE,
    NDE = NDE, 
    TE = TE,
    PTE = PTE
  ))
}


# Function to generate survival and SLD data
generate_data <- function(params, n_patients, link_function = "current_sld", seed_num = 4) {
  Z <- generate_block_randomization(n_patients)
  time_points_per_patient <- lapply(1:n_patients, function(x) generate_time_points())
  
  # SLD data
  sld_data <- simulate_sld(
    n = n_patients,
    mu_M0 = params$mu_M0,
    mu_kg = params$mu_kg,
    mu_ks = params$mu_ks,
    omega_M0 = params$omega_M0,
    omega_kg = params$omega_kg,
    omega_ks = params$omega_ks,
    beta = params$beta,
    Z = Z,
    time_points_per_patient = time_points_per_patient,
    seed_num = seed_num
  )
  
  # Apply 20% rule
  sld_data_filtered <- apply_sld_rule(sld_data, seed_num = seed_num)
  
  # Survival data
  survival_results <- simulate_survival_jm_weibull(
    n = n_patients,
    lambda = params$lambda,
    kappa = params$kappa,
    gamma = params$gamma,
    beta_os = params$beta_os,
    Z = Z,
    time_points_per_patient = lapply(1:n_patients, function(x) seq(from = 0, to = 2, by = .001)),
    ks_params = sld_data |> distinct(id, ks_i) |> pull(ks_i),
    kg_params = sld_data |> distinct(id, kg_i) |> pull(kg_i),
    M0_values = sld_data |> distinct(id, M0_i) |> pull(M0_i),
    link_function = link_function,
    beta = params$beta,
    seed_num = seed_num
  )
  
  # Combine the survival times with the filtered SLD dataset
  data_os <- data.frame(
    id = 1:n_patients,
    Z = Z,
    survival_times = survival_results$survival_times,
    status = survival_results$event_observed
  )
  
  sld_data_filtered <- sld_data_filtered %>%
    left_join(data_os, by = c("id", "Z")) %>%
    filter(time <= survival_times) %>%
    dplyr::select(-survival_times, -status) # Remove columns after filtering
  
  
  list(sld = sld_data_filtered, survival = data_os)
}

# Function to generate time points for each patient
generate_time_points <- function() {
  time_points_6_weeks <- seq(0, 365, by = 42)
  time_points_9_weeks <- seq(365 + 63, 2 * 365, by = 63)
  time_points <- c(time_points_6_weeks, time_points_9_weeks)
  time_points <- time_points + sample(0:8, length(time_points), replace = TRUE)
  time_points[1] <- sample(-14:-1, 1)
  sort(time_points / 365)
}



# Function to prepare Stan data
prepare_stan_data <- function(sld_data_filtered, survival_results, params) {
  # Prepare data similar to your existing code
  data_os <- survival_results
  tumor.tot <- sld_data_filtered |>
    mutate(
      USUBJID = paste0("subject_", id),
      TIME = time,
      SLD = sld,
      ARM.PRED = factor(Z, levels = c(0, 1), labels = c( "0", "1"), ordered = T)
    )
  
  ate <- data_os |>
    mutate(
      AYR = survival_times,
      USUBJID = paste0("subject_", id),
      CNSR = 1 - status
    )
  
  pts_ls <- tumor.tot |> distinct(USUBJID)
  
  data_sld <- tumor.tot |>
    ungroup() |>
    select(-id) |>
    rename(TRT01P = ARM.PRED) |>
    mutate(
      AYR = TIME,
      all_arms = "all_arms",
      SLDC = SLD
    ) |>
    filter(USUBJID %in% pts_ls$USUBJID)
  
  # Create column with the time of the last observed measurement
  last_obs_time_df <- data_sld |>
    ungroup() |>
    group_by(USUBJID) |>
    arrange(USUBJID, -AYR) |>
    dplyr::select(USUBJID, SLDC, AYR) |>
    slice(1) |>
    rename(
      last_obs_time = AYR,
      last_sld_obs = SLDC
    )
  
  data_sld <- data_sld |>
    left_join(last_obs_time_df, by = "USUBJID") |>
    mutate(id = as.numeric(factor(USUBJID))) |>
    arrange(id, AYR)
  
  
  
  
  data_os <- ate |>
    select(-id) |>
    inner_join(data_sld[, c("USUBJID", "TRT01P", "id")], by = "USUBJID") |>
    distinct(USUBJID, .keep_all = T) |>
    filter(USUBJID %in% data_sld$USUBJID) |>
    arrange(id)
  
  ## Extract sparse parts
  mat_inds_obs_y_all <- t(model.matrix(~ -1 + USUBJID, data = data_sld))
  sparse_mat_inds_obs_y_all <- extract_sparse_parts(mat_inds_obs_y_all)
  
  # Censoring threshold for tumor size
  cens_threshold <- max(data_sld |> select(USUBJID, SLD) |> group_by(USUBJID) |> arrange(SLD) |> slice(1) |> arrange(SLD) |> ungroup()|> slice(2) |> pull(), 2.5)
  
  # Derive sparse matrices
  obs_y_dat <- subset(data_sld, SLDC > cens_threshold)
  mat_inds_obs_y <- t(model.matrix(~ -1 + USUBJID, data = obs_y_dat))
  sparse_mat_inds_obs_y <- extract_sparse_parts(mat_inds_obs_y)
  
  cens_y_dat <- subset(data_sld, SLDC <= cens_threshold)
  mat_inds_cens_y <- t(model.matrix(~ -1 + USUBJID, data = cens_y_dat))
  sparse_mat_inds_cens_y <- extract_sparse_parts(mat_inds_cens_y)
  
  # Timepoints for OS hazard and survival function estimation to generate predictions
  os_pred_times <- c(.1, .5, 1, 1.5, 2)
  sld_exp_times <- c(.1, .5, 1, 1.5, 2)
  
  # 0 should be control
  dose_per_pts <- data_sld |>
    dplyr::select(USUBJID, TRT01P) |>
    distinct(USUBJID, .keep_all = T) |>
    dplyr::select(TRT01P) |>
    mutate(TRT01P = as.numeric(factor(TRT01P, levels = c("0", "1"))) - 1)
  
  
  gh_parameters <- statmod::gauss.quad(n = 15, kind = "legendre")
  
  ## Overall survival covariates design matrix
  os_cov_design <- model.matrix(~ TRT01P - 1, data = data_os)[, 2, drop = FALSE]
  
  s0_design <- os_cov_design
  s0_design <- matrix(0, ncol = ncol(os_cov_design), nrow = nrow(os_cov_design))
  s1_design <- os_cov_design
  s1_design <- matrix(1, ncol = ncol(os_cov_design), nrow = nrow(os_cov_design))
  
  stan_data <- list(
    Nind = length(unique(data_sld$USUBJID)), #  Number of individuals.
    Nind_cens = nrow(distinct(cens_y_dat, USUBJID)),
    Nta_total = nrow(drop_na(data_sld, SLDC)), # Total number of tumor assessments.
    Nta_obs_y = sum(data_sld$SLD > cens_threshold),
    ind_index = data_sld$id, # Index of individuals for each tumor assessment.
    obs_y_index = which(data_sld$SLD > cens_threshold), # Index of tumor assessments .
    cens_y_index = which(data_sld$SLD <= cens_threshold), # Index of samll tumor assessments .
    
    Nta_cens_y = sum(data_sld$SLD <= cens_threshold),
    Yobs = data_sld$SLDC, # Array of individual responses.
    Tobs = data_sld$AYR,
    Ythreshold = cens_threshold,
    Tobs_neg = data_sld |> distinct(USUBJID, .keep_all = T) |> group_by(USUBJID) |>
      arrange(AYR) |> slice(1) |> ungroup() |> arrange(id) |> dplyr::select(AYR) |> unlist(),
    Nind_dead = sum(data_os$CNSR == 0),
    CNSR = data_os$CNSR,
    
    ## Treatment related data
    trt = as.numeric(factor(data_sld$TRT01P, levels = c("0", "1"))) - 1,
    trt_gen_q = data_sld |> ungroup() |> distinct(USUBJID, .keep_all = T) |> arrange(id) |> dplyr::select(TRT01P) |>
      unlist() |> factor(levels = c("0", "1")) |> as.numeric() - 1,
    trt_pts = dose_per_pts$TRT01P,
    
    # Matrix of individuals x observed tumor assessments.
    n_w_mat_inds_obs_y_all = length(sparse_mat_inds_obs_y_all$w),
    w_mat_inds_obs_y_all = sparse_mat_inds_obs_y_all$w,
    n_v_mat_inds_obs_y_all = length(sparse_mat_inds_obs_y_all$v),
    v_mat_inds_obs_y_all = sparse_mat_inds_obs_y_all$v,
    n_u_mat_inds_obs_y_all = length(sparse_mat_inds_obs_y_all$u),
    u_mat_inds_obs_y_all = sparse_mat_inds_obs_y_all$u,
    
    
    # Matrix of individuals x observed tumor assessments.
    n_w_mat_inds_obs_y = length(sparse_mat_inds_obs_y$w),
    w_mat_inds_obs_y = sparse_mat_inds_obs_y$w,
    n_v_mat_inds_obs_y = length(sparse_mat_inds_obs_y$v),
    v_mat_inds_obs_y = sparse_mat_inds_obs_y$v,
    n_u_mat_inds_obs_y = length(sparse_mat_inds_obs_y$u),
    u_mat_inds_obs_y = sparse_mat_inds_obs_y$u,
    
    # Matrix of individuals x censored tumor assessments.
    n_w_mat_inds_cens_y = length(sparse_mat_inds_cens_y$w),
    w_mat_inds_cens_y = sparse_mat_inds_cens_y$w,
    n_v_mat_inds_cens_y = length(sparse_mat_inds_cens_y$v),
    v_mat_inds_cens_y = sparse_mat_inds_cens_y$v,
    n_u_mat_inds_cens_y = length(sparse_mat_inds_cens_y$u),
    u_mat_inds_cens_y = sparse_mat_inds_cens_y$u,
    sld_exp_times = sld_exp_times,
    sld_pred_times = sld_exp_times, # must be equal to sld_exp_times
    final_sld_time = unlist(distinct(data_sld, USUBJID, .keep_all = T) |>
                              arrange(id) |>
                              dplyr::select(last_obs_time)),
    n_sld_exp_times = length(sld_exp_times),
    n_sld_pred_times = length(sld_exp_times),
    final_sld = unlist(distinct(data_sld, USUBJID, .keep_all = T) |>
                         dplyr::select(last_sld_obs)),
    total_obs_sld_times = length(data_sld$SLDC),
    n_arms = data_sld |> distinct(all_arms) |> unlist() |> length(),
    
    # The patients in each of the four different treatment arms.
    n_index_per_arm = c(
      data_sld |>
        distinct(USUBJID, .keep_all = T) |>
        group_by(TRT01P) |>
        tally() |>
        dplyr::select(n) |>
        unlist()
    ),
    index_per_arm = as.numeric(factor(data_sld |>
                                        distinct(USUBJID, .keep_all = T) |> arrange(id) |>
                                        dplyr::select(all_arms) |> unlist())),
    arm_index = as.numeric(factor(data_sld |>
                                    distinct(USUBJID, .keep_all = T) |> arrange(id) |>
                                    dplyr::select(all_arms) |> unlist())),
    arm_index_real = as.numeric(factor(data_sld |>
                                         distinct(USUBJID, .keep_all = T) |> arrange(id) |>
                                         dplyr::select(TRT01P) |> unlist())),
    
    # Survival data.
    Times = data_os$AYR,
    os_ind_index = data_os$id,
    os_cov_design = os_cov_design,
    s0_design = s0_design,
    s1_design = s1_design,
    p_os_cov_design = ncol(os_cov_design),
    dead_ind_index = data_os[data_os$CNSR == 0, ]$id,
    n_nodes = length(gh_parameters$nodes),
    nodes = gh_parameters$nodes,
    weights = gh_parameters$weights,
    
    # For OS predictions.
    os_pred_times = os_pred_times,
    n_os_pred_times = length(os_pred_times)
  )
  return(stan_data)
}


# Function to calculate survival probabilities based on ks and Z
calculate_effects <- function(time_points, lambda, kappa, gamma, beta_os, link_function, current_k_s, k_g, M0, current_Z) {
  sapply(time_points, function(t) {
    # Link function based on time and ks
    link_value <- link_function(t, current_k_s, k_g, M0)
    # Cumulative hazard up to time t
    H_t <- (t / lambda)^kappa
    H_t <- H_t * exp(gamma * link_value + beta_os * current_Z)
    # Survival probability as exponent of negative cumulative hazard
    exp(-H_t)
  })
}


generate_inits <- function(stan_data, mu_bsld_init = 60, mu_ks_init = 0.1, mu_kg_init = 0.1, 
                           omega_init = 0.1, beta_trt_init = 0, sigma_init = 1, 
                           lambda_init = 1, p_init = 2, gamma_init = 0, beta_os_cov_init = 0) {
  
  # Create a list of initial values based on the input parameters and stan_data structure
  inits <- list(
    list(
      "mu_bsld" = mu_bsld_init,
      "mu_ks" = mu_ks_init,
      "mu_kg" = mu_kg_init,
      "omega_bsld" = omega_init,
      "omega_ks" = omega_init,
      "omega_kg" = omega_init,
      "eta_tilde_bsld" = rep(0, stan_data$Nind),
      "eta_tilde_ks" = rep(0, stan_data$Nind),
      "eta_tilde_kg" = rep(0, stan_data$Nind),
      "mean_mu_ks" = rep(0.01, stan_data$n_arms),
      "mean_mu_kg" = rep(0.1, stan_data$n_arms),
      "mean_mu_bsld" = rep(10, stan_data$n_arms),
      "sd_mu_ks" = rep(1, stan_data$n_arms),
      "sd_mu_kg" = rep(1, stan_data$n_arms),
      "sd_mu_bsld" = rep(0.01, stan_data$n_arms),
      "beta_trt" = beta_trt_init,
      "sigma" = sigma_init,
      "lambda" = lambda_init,
      "p" = p_init,
      "gamma" = gamma_init,
      "beta_os_cov" = rep(beta_os_cov_init, ncol(stan_data$os_cov_design))
    ),
    list(
      "mu_bsld" = mu_bsld_init * 0.8,
      "mu_ks" = mu_ks_init * 0.1,
      "mu_kg" = mu_kg_init * 1.1,
      "omega_bsld" = omega_init * 1.1,
      "omega_ks" = omega_init * 1.1,
      "omega_kg" = omega_init * 1.1,
      "eta_tilde_bsld" = rep(0, stan_data$Nind),
      "eta_tilde_ks" = rep(0, stan_data$Nind),
      "eta_tilde_kg" = rep(0, stan_data$Nind),
      "mean_mu_ks" = rep(0.01, stan_data$n_arms),
      "mean_mu_kg" = rep(0.1, stan_data$n_arms),
      "mean_mu_bsld" = rep(10, stan_data$n_arms),
      "sd_mu_ks" = rep(1, stan_data$n_arms),
      "sd_mu_kg" = rep(1, stan_data$n_arms),
      "sd_mu_bsld" = rep(0.01, stan_data$n_arms),
      "beta_trt" = beta_trt_init * 0.5,
      "sigma" = sigma_init * 1.1,
      "lambda" = lambda_init * 1.1,
      "p" = p_init - 0.5,
      "gamma" = gamma_init - 0.1,
      "beta_os_cov" = rep(beta_os_cov_init * 1.5, ncol(stan_data$os_cov_design))
    ),
    list(
      "mu_bsld" = mu_bsld_init * 1.1,
      "mu_ks" = mu_ks_init * 0.9,
      "mu_kg" = mu_kg_init * 1.2,
      "omega_bsld" = omega_init * 1.2,
      "omega_ks" = omega_init * 1.2,
      "omega_kg" = omega_init * 1.2,
      "eta_tilde_bsld" = rep(0.5, stan_data$Nind),
      "eta_tilde_ks" = rep(0.5, stan_data$Nind),
      "eta_tilde_kg" = rep(0.5, stan_data$Nind),
      "mean_mu_ks" = rep(0.02, stan_data$n_arms),
      "mean_mu_kg" = rep(0.2, stan_data$n_arms),
      "mean_mu_bsld" = rep(15, stan_data$n_arms),
      "sd_mu_ks" = rep(1.5, stan_data$n_arms),
      "sd_mu_kg" = rep(1.5, stan_data$n_arms),
      "sd_mu_bsld" = rep(0.02, stan_data$n_arms),
      "beta_trt" = beta_trt_init * -0.5,
      "sigma" = sigma_init * 2,
      "lambda" = lambda_init * 0.9,
      "p" = p_init * 1.1,
      "gamma" = gamma_init * 1.5,
      "beta_os_cov" = rep(beta_os_cov_init * 0.8, ncol(stan_data$os_cov_design))
    )
  )
  
  return(inits)
}

generate_inits_ttn <- function(stan_data, mu_bsld_init = 60, mu_ks_init = 0.1, mu_kg_init = 0.1, 
                           omega_init = 0.1, beta_trt_init = 0, sigma_init = 1, 
                           lambda_init = 1, p_init = 2, gamma_init = 0, beta_os_cov_init = 0) {
  
  # Create a list of initial values based on the input parameters and stan_data structure
  inits <- list(
    list(
      "mu_bsld" = mu_bsld_init,
      "mu_ks" = mu_ks_init,
      "mu_kg" = mu_kg_init,
      "omega_bsld" = omega_init,
      "omega_ks" = omega_init,
      "omega_kg" = omega_init,
      "eta_tilde_bsld" = rep(0, stan_data$Nind),
      "eta_tilde_ks" = rep(0, stan_data$Nind),
      "eta_tilde_kg" = rep(0, stan_data$Nind),
      "mean_mu_ks" = rep(0.01, stan_data$n_arms),
      "mean_mu_kg" = rep(0.1, stan_data$n_arms),
      "mean_mu_bsld" = rep(10, stan_data$n_arms),
      "sd_mu_ks" = rep(1, stan_data$n_arms),
      "sd_mu_kg" = rep(1, stan_data$n_arms),
      "sd_mu_bsld" = rep(0.01, stan_data$n_arms),
      "beta_trt" = beta_trt_init,
      "sigma" = sigma_init,
      "lambda" = lambda_init,
      "p" = p_init,
      "gamma" = gamma_init,
      "beta_os_cov" = rep(beta_os_cov_init, ncol(stan_data$os_cov_design)),
      "tau" = 1,
      "lambda_os" = rep(.3, ncol(stan_data$os_cov_design)),
      "lambda_trt" = 1,
      "lambda_gamma" = .82,
      "z_os" = rep(-.1, ncol(stan_data$os_cov_design)),
      "z_trt" = .3,
      "z_gamma" = -.2
    ),
    list(
      "mu_bsld" = mu_bsld_init * 0.8,
      "mu_ks" = mu_ks_init * 0.1,
      "mu_kg" = mu_kg_init * 1.1,
      "omega_bsld" = omega_init * 1.1,
      "omega_ks" = omega_init * 1.1,
      "omega_kg" = omega_init * 1.1,
      "eta_tilde_bsld" = rep(0, stan_data$Nind),
      "eta_tilde_ks" = rep(0, stan_data$Nind),
      "eta_tilde_kg" = rep(0, stan_data$Nind),
      "mean_mu_ks" = rep(0.01, stan_data$n_arms),
      "mean_mu_kg" = rep(0.1, stan_data$n_arms),
      "mean_mu_bsld" = rep(10, stan_data$n_arms),
      "sd_mu_ks" = rep(1, stan_data$n_arms),
      "sd_mu_kg" = rep(1, stan_data$n_arms),
      "sd_mu_bsld" = rep(0.01, stan_data$n_arms),
      "beta_trt" = beta_trt_init * 0.5,
      "sigma" = sigma_init * 1.1,
      "lambda" = lambda_init * 1.1,
      "p" = p_init - 0.5,
      "gamma" = gamma_init - 0.1,
      "beta_os_cov" = rep(beta_os_cov_init * 1.5, ncol(stan_data$os_cov_design)),
      "tau" = 0.8,
      "lambda_os" = rep(0.5, ncol(stan_data$os_cov_design)),
      "lambda_trt" = .7,
      "lambda_gamma" = .3,
      "z_os" = rep(-.5, ncol(stan_data$os_cov_design)),
      "z_trt" = 0.4,
      "z_gamma" = -.5
    ),
    list(
      "mu_bsld" = mu_bsld_init * 1.1,
      "mu_ks" = mu_ks_init * 0.9,
      "mu_kg" = mu_kg_init * 1.2,
      "omega_bsld" = omega_init * 1.2,
      "omega_ks" = omega_init * 1.2,
      "omega_kg" = omega_init * 1.2,
      "eta_tilde_bsld" = rep(0.5, stan_data$Nind),
      "eta_tilde_ks" = rep(0.5, stan_data$Nind),
      "eta_tilde_kg" = rep(0.5, stan_data$Nind),
      "mean_mu_ks" = rep(0.02, stan_data$n_arms),
      "mean_mu_kg" = rep(0.2, stan_data$n_arms),
      "mean_mu_bsld" = rep(15, stan_data$n_arms),
      "sd_mu_ks" = rep(1.5, stan_data$n_arms),
      "sd_mu_kg" = rep(1.5, stan_data$n_arms),
      "sd_mu_bsld" = rep(0.02, stan_data$n_arms),
      "beta_trt" = beta_trt_init * -0.5,
      "sigma" = sigma_init * 2,
      "lambda" = lambda_init * 0.9,
      "p" = p_init * 1.1,
      "gamma" = gamma_init * 1.5,
      "beta_os_cov" = rep(beta_os_cov_init * 0.8, ncol(stan_data$os_cov_design)),
      "tau" = 0.5,
      "lambda_os" = rep(.5, ncol(stan_data$os_cov_design)),
      "lambda_trt" = .5,
      "lambda_gamma" = 1,
      "z_os" = rep(-.7, ncol(stan_data$os_cov_design)),
      "z_trt" = 0.5,
      "z_gamma" = -1
    )
  )
  
  return(inits)
}


get_link_value <- function(t, k_s, k_g, M0, link_function = "current_sld") {
  if (link_function == "current_sld") {
    if (t < 0) {
      return(M0 * exp(k_g * t))  # Before treatment: only growth
    } else {
      return(M0 * (exp(k_g * t) + exp(-k_s * t) - 1))  # After treatment: growth + shrinkage
    }
    
  } else if (link_function == "sld_slope") {
    return(M0 * (k_g * exp(k_g * t) - k_s * exp(-k_s * t)))
    
  } else if (link_function == "ks_parameter") {
    return(log(k_s))
    
  } else if (link_function == "time_to_nadir") {
    t_nadir <- log(k_s / k_g) / (k_g + k_s)
    return(t_nadir)
    
  } else {
    stop("Unknown link function specified.")
  }
}
