# List of required packages
required_packages <- c(
  "dplyr", "truncnorm", "survival", "cmdstanr", 
  "posterior", "ggplot2", "ggpubr"
)

# Install any missing packages
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg)
  }
}

# Load Libraries
library(dplyr)
library(truncnorm)
library(survival)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(ggpubr)
cmdstanr::check_cmdstan_toolchain()
install_cmdstan()

# Load custom functions
source("Functions.R")

# Compile Stan model
stan_model <- cmdstan_model(stan_file = "gq_tgi_beta_ks_os_trt_surv_weibull_link_ks.stan")


# Define parameters for each scenario
scenarios <- list(
  scenario_1 = list(mu_M0 = 60, mu_kg = 0.35, mu_ks = 0.53, beta = 0.5, lambda = 1.2, kappa = 1.5, gamma = -0.2, beta_os = -0.5, omega_ks = 1.3, omega_kg = 1, omega_M0 = .6),
  scenario_2 = list(mu_M0 = 60, mu_kg = 0.35, mu_ks = 0.53, beta = .7, lambda = 1.2, kappa = 1.5, gamma = -1, beta_os = -0.05, omega_ks = 1.3, omega_kg = 1, omega_M0 = .6),
  scenario_3 = list(mu_M0 = 60, mu_kg = 0.35, mu_ks = 0.53, beta = 0.1, lambda = 1.2, kappa = 1.5, gamma = -0.05, beta_os = -0.7, omega_ks = 1.3, omega_kg = 1, omega_M0 = .6)
)

# List of sample sizes for small and large datasets
sample_sizes <- list(small = 200, large = 800)

# Define parameters to extract
params_to_extract <- c("beta_os_cov", "gamma", "p", "lambda", "beta_trt",
                       "mu_kg[1]", "mu_ks[1]", "mu_bsld", "sigma", 
                       "omega_kg", "omega_ks", "omega_bsld")
n_reps <- 100

for (rep in 1:n_reps) {
  cat("Replication:", rep, "\n")
  
  # Store datasets and results for each replication
  replication_results <- list()
  inits_list <- list()
  for (scenario_name in names(scenarios)) {
    for (size in names(sample_sizes)) { 
      params <- scenarios[[scenario_name]]
      n_patients <- sample_sizes[[size]]
      
      # Generate SLD and survival data
      generated_data <- generate_data(params, n_patients, link_function = "ks_parameter", seed_num = rep)
      sld_data_filtered <- generated_data$sld
      survival_data <- generated_data$survival
      
      # Prepare Stan data
      stan_data <- prepare_stan_data(sld_data_filtered, survival_data, params)
      inits_list[[paste0("stan_data_", scenario_name, "_", size)]] <- generate_inits(stan_data)
      # Fit the Stan model
      stan_fit <- stan_model$sample(
        data = stan_data,
        iter_sampling = 1000,
        iter_warmup = 2000,
        parallel_chains = 3,
        chains = 3,
        init = inits_list[[paste0("stan_data_", scenario_name, "_", size)]]
      )
      
      
      # Extract posterior draws for the parameters
      posterior_draws <- stan_fit$draws(variables = params_to_extract, format = "matrix")
      
      # Calculate 95% credible intervals
      credible_intervals <- apply(posterior_draws, 2, function(x) {
        quantile(x, probs = c(0.025, 0.975))
      })
      
      # Convert to a data frame for better readability
      credible_intervals_df <- as.data.frame(t(credible_intervals))
      colnames(credible_intervals_df) <- c("Lower_95", "Upper_95")
      
      # Example of combining with other summary statistics
      summary_stats <- stan_fit$summary(variables = params_to_extract)
      summary_stats <- cbind(summary_stats, credible_intervals_df)
      summary_medians <- summary_stats 
      
      
      # Extract log-likelihood matrix
      log_lik <- stan_fit$draws(variables = "log_lik", format = "matrix")
      
      # Convert to LOO object
      loo_object <- loo::loo(log_lik)
      
      # Save LOOIC
      looic_value <- loo_object$estimates["looic", "Estimate"]
      
      posterior_samples <- stan_fit$draws(variables = c("NIE_pop_samples", "NDE_pop_samples", "PTE_pop_samples"))
      posterior_df <- as_draws_df(posterior_samples)
      
      
      posterior_long <- posterior_df %>%
        pivot_longer(
          cols = starts_with("NIE_pop_samples") | starts_with("NDE_pop_samples") | starts_with("PTE_pop_samples"),
          names_to = c("Metric", "Time"),
          names_pattern = "(.*)\\[(\\d+)\\]",
          values_to = "Value"
        ) %>%
        mutate(Time = as.numeric(Time))
      
      effects_summary_stats <- posterior_long %>%
        group_by(Metric, Time) %>%
        summarize(
          Median = median(Value),
          SD = sd(Value),
          Lower_95 =  quantile(Value, probs = c(0.025), na.rm = TRUE),
          Upper_95 =  quantile(Value, probs = c(0.975), na.rm = TRUE),
          .groups = "drop"
        )
      
      
      posterior_samples <- stan_fit$draws(variables = c("rizopoulos_marginal_NIE", "rizopoulos_marginal_NDE", "rizopoulos_marginal_PTE"))
      posterior_df <- as_draws_df(posterior_samples)
      
      
      posterior_long <- posterior_df %>%
        pivot_longer(
          cols = starts_with("rizopoulos_marginal_NIE") | starts_with("rizopoulos_marginal_NDE") | starts_with("rizopoulos_marginal_PTE"),
          names_to = c("Metric", "Time"),
          names_pattern = "(.*)\\[(\\d+)\\]",
          values_to = "Value"
        ) %>%
        mutate(Time = as.numeric(Time))
      
      rizopoulos_effects_summary_stats <- posterior_long %>%
        group_by(Metric, Time) %>%
        summarize(
          Median = median(Value),
          SD = sd(Value),
          Lower_95 =  quantile(Value, probs = c(0.025), na.rm = T),
          Upper_95 =  quantile(Value, probs = c(0.975), na.rm = T),
          .groups = "drop"
        )
      
      
      posterior_samples_typical <- stan_fit$draws(variables = c("typical_NIE", "typical_NDE", "typical_PTE"))
      posterior_df_typical <- as_draws_df(posterior_samples_typical)
      
      posterior_long_typical <- posterior_df_typical %>%
        pivot_longer(
          cols = starts_with("typical_NIE") | starts_with("typical_NDE") | starts_with("typical_PTE"),
          names_to = c("Metric", "Time"),
          names_pattern = "(.*)\\[(\\d+)\\]",
          values_to = "Value"
        ) %>%
        mutate(Time = as.numeric(Time))
      
      effects_typical_summary_stats <- posterior_long_typical %>%
        group_by(Metric, Time) %>%
        summarize(
          Median = median(Value),
          SD = sd(Value),
          Lower_95 =  quantile(Value, probs = c(0.025)),
          Upper_95 =  quantile(Value, probs = c(0.975)),
          .groups = "drop"
        )
      
      
      # Combine medians and SEs into a single data frame for the current model
      summary_combined <- summary_medians 
      
      # Store results for this scenario and sample size
      replication_results[[paste0(scenario_name, "_", size)]] <- list(
        stan_fit = stan_fit,
        sld_data = sld_data_filtered,
        survival_data = survival_data,
        summary = summary_combined,
        looic = looic_value,
        effects = effects_summary_stats,
        rizopoulos_effects_marginal = rizopoulos_effects_summary_stats,
        effects_typical = effects_typical_summary_stats
      )
      
    }
  }
  
  # Save each replication's results
  all_results[[paste0("rep_", rep)]] <- replication_results
}

# Placeholder to store the results
summary_results <- list()

# Define custom rounding function with selective precision
round_values <- function(df) {
  df %>%
    mutate(
      True_Value = round(True_Value, 4),
      Mean_of_Medians = round(Mean_of_Medians, 4),
      Mean_Bias = round(Mean_Bias, 3),
      SD_Bias = round(SD_Bias, 4),
      Relative_Bias_Percent = round(Relative_Bias_Percent, 1),
      Coverage_Percentage = round(Coverage_Percentage, 1),
      SD_of_Medians = round(SD_of_Medians, 3),
      Mean_of_SDs = round(Mean_of_SDs, 3)
    )
}

# Time points for PTE, NIE, and NDE
time_points <- c(seq(from = 0, to = 2, by = 0.1))


############ True Marginal effects ##################


# Time points for PTE, NIE, and NDE
time_points <- c(seq(from = 0, to = 2, by = 0.1))

# Placeholder for true effects across all replications
true_marginal_effects_results <- list()

# Loop over each scenario and sample size
for (scenario_name in names(scenarios)) {
  for (size in names(sample_sizes)) {
    # Placeholder to store results for this scenario and size
    scenario_effects <- data.frame(time = numeric(), NIE = numeric(), NDE = numeric(), PTE = numeric(), replication = integer())
    
    # Loop through each replication
    for (rep in 1:n_reps) {
      # Extract the generated dataset for this replication
      generated_data <- all_results[[paste0("rep_", rep)]][[paste0(scenario_name, "_", size)]]
      sld_data <- generated_data$sld_data
      
      # Extract random effects per individual
      M0_i <- scenarios[[scenario_name]]$mu_M0 * exp(rnorm(10000, 0,scenarios[[scenario_name]]$omega_M0))
      k_g_i <- scenarios[[scenario_name]]$mu_kg * exp(rnorm(10000, 0,scenarios[[scenario_name]]$omega_kg))
      k_s_i <- scenarios[[scenario_name]]$mu_ks * exp(rnorm(10000, 0, scenarios[[scenario_name]]$omega_ks))
      
      # Initialize vectors to store results for this replication
      NIE_rep <- numeric(length(time_points))
      NDE_rep <- numeric(length(time_points))
      PTE_rep <- numeric(length(time_points))
      
      # Loop over each time point
      for (t_idx in seq_along(time_points)) {
        t <- time_points[t_idx]
        
        # Calculate survival probabilities for each individual
        S_00 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i, k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 0))
        S_01 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i, k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 1))
        S_10 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i * exp(scenarios[[scenario_name]]$beta), k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 0))
        S_11 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i * exp(scenarios[[scenario_name]]$beta), k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 1))
        
        # Calculate NIE, NDE, and PTE for each individual
        NIE_ind <- S_11 - S_01
        NDE_ind <- S_01 - S_00
        PTE_ind <- ifelse((S_11 - S_00) != 0, NIE_ind / (S_11 - S_00), NA)
        
        # Aggregate across individuals (e.g., take the median)
        NIE_rep[t_idx] <- median(NIE_ind, na.rm = TRUE)
        NDE_rep[t_idx] <- median(NDE_ind, na.rm = TRUE)
        PTE_rep[t_idx] <- median(PTE_ind, na.rm = TRUE)
      }
      
      # Store results for this replication
      scenario_effects <- rbind(scenario_effects, data.frame(
        time = time_points,
        NIE = NIE_rep,
        NDE = NDE_rep,
        PTE = PTE_rep,
        replication = rep
      ))
    }
    
    # Aggregate across replications
    final_effects <- scenario_effects %>%
      group_by(time) %>%
      summarize(
        Median_NIE = median(NIE, na.rm = TRUE),
        Median_NDE = median(NDE, na.rm = TRUE),
        Median_PTE = median(PTE, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Save results for this scenario and size
    true_marginal_effects_results[[paste0(scenario_name, "_", size)]] <- final_effects
  }
}

# Combine results into a single data frame for all scenarios and sizes
combined_true_marginal_effects <- bind_rows(true_marginal_effects_results, .id = "Scenario_SampleSize") |> 
  rename(NIE = Median_NIE, NDE = Median_NDE, PTE = Median_PTE ) |> 
  mutate(Scenario_Sample_Size = Scenario_SampleSize)

combined_true_marginal_effects$Scenario_Sample_Size <- factor(
  combined_true_marginal_effects$Scenario_Sample_Size,
  levels = c(
    "scenario_1_small", 
    "scenario_2_small", 
    "scenario_3_small", 
    "scenario_1_large", 
    "scenario_2_large", 
    "scenario_3_large"
  )
)

############ True Typical effects ##################


# Time points for PTE, NIE, and NDE
time_points <- c(seq(from = 0, to = 2, by = 0.1))

# Placeholder for true effects across all replications
true_typical_effects_results <- list()

# Loop over each scenario and sample size
for (scenario_name in names(scenarios)) {
  for (size in names(sample_sizes)) {
    # Placeholder to store results for this scenario and size
    scenario_effects <- data.frame(time = numeric(), NIE = numeric(), NDE = numeric(), PTE = numeric(), replication = integer())
    
    # Loop through each replication
    for (rep in 1:n_reps) {
      # Extract the generated dataset for this replication
      generated_data <- all_results[[paste0("rep_", rep)]][[paste0(scenario_name, "_", size)]]
      sld_data <- generated_data$sld_data
      
      # Extract random effects per individual
      M0_i <- scenarios[[scenario_name]]$mu_M0
      k_g_i <- scenarios[[scenario_name]]$mu_kg
      k_s_i <- scenarios[[scenario_name]]$mu_ks
      
      # Initialize vectors to store results for this replication
      NIE_rep <- numeric(length(time_points))
      NDE_rep <- numeric(length(time_points))
      PTE_rep <- numeric(length(time_points))
      
      # Loop over each time point
      for (t_idx in seq_along(time_points)) {
        t <- time_points[t_idx]
        
        # Calculate survival probabilities for each individual
        S_00 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i, k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 0))
        S_01 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i, k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 1))
        S_10 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i * exp(scenarios[[scenario_name]]$beta), k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 0))
        S_11 <- exp(-((t / scenarios[[scenario_name]]$lambda)^scenarios[[scenario_name]]$kappa) *
                      exp(scenarios[[scenario_name]]$gamma * get_link_value(t, k_s_i * exp(scenarios[[scenario_name]]$beta), k_g_i, M0_i, link_function = "ks_parameter") + scenarios[[scenario_name]]$beta_os * 1))
        
        # Calculate NIE, NDE, and PTE for each individual
        NIE_ind <- S_11 - S_01
        NDE_ind <- S_01 - S_00
        PTE_ind <- ifelse((S_11 - S_00) != 0, NIE_ind / (S_11 - S_00), NA)
        
        # Aggregate across individuals (e.g., take the median)
        NIE_rep[t_idx] <- median(NIE_ind, na.rm = TRUE)
        NDE_rep[t_idx] <- median(NDE_ind, na.rm = TRUE)
        PTE_rep[t_idx] <- median(PTE_ind, na.rm = TRUE)
      }
      
      # Store results for this replication
      scenario_effects <- rbind(scenario_effects, data.frame(
        time = time_points,
        NIE = NIE_rep,
        NDE = NDE_rep,
        PTE = PTE_rep,
        replication = rep
      ))
    }
    
    # Aggregate across replications
    final_effects <- scenario_effects %>%
      group_by(time) %>%
      summarize(
        Median_NIE = median(NIE, na.rm = TRUE),
        Median_NDE = median(NDE, na.rm = TRUE),
        Median_PTE = median(PTE, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Save results for this scenario and size
    true_typical_effects_results[[paste0(scenario_name, "_", size)]] <- final_effects
  }
}

# Combine results into a single data frame for all scenarios and sizes
combined_true_typical_effects <- bind_rows(true_typical_effects_results, .id = "Scenario_SampleSize") |> 
  rename(NIE = Median_NIE, NDE = Median_NDE, PTE = Median_PTE ) |> 
  mutate(Scenario_Sample_Size = Scenario_SampleSize)

combined_true_typical_effects$Scenario_Sample_Size <- factor(
  combined_true_typical_effects$Scenario_Sample_Size,
  levels = c(
    "scenario_1_small", 
    "scenario_2_small", 
    "scenario_3_small", 
    "scenario_1_large", 
    "scenario_2_large", 
    "scenario_3_large"
  )
)



################## Summary stats ##################


time_points <- c(0.5, 1, 2)
# Loop over each scenario and sample size
for (scenario_name in names(scenarios)) {
  for (size in names(sample_sizes)[1]) {
    # Initialize data frames for medians, biases, and coverage
    medians_df <- data.frame(matrix(NA, nrow = n_reps, ncol = length(params_to_extract)))
    colnames(medians_df) <- params_to_extract
    SDs_df <- data.frame(matrix(NA, nrow = n_reps, ncol = length(params_to_extract)))
    colnames(SDs_df) <- params_to_extract
    bias_df <- data.frame(matrix(NA, nrow = n_reps, ncol = length(params_to_extract)))
    colnames(bias_df) <- params_to_extract
    rhat_df <- data.frame(matrix(NA, nrow = n_reps, ncol = length(params_to_extract)))
    colnames(rhat_df) <- params_to_extract
    coverage_flags <- data.frame(matrix(NA, nrow = n_reps, ncol = length(params_to_extract)))
    colnames(coverage_flags) <- params_to_extract
    rhat_flags <- data.frame(matrix(NA, nrow = n_reps, ncol = length(params_to_extract)))
    colnames(rhat_flags) <- params_to_extract
    coverage_flags_effects <- data.frame(matrix(NA, nrow = n_reps, ncol = 3*length(time_points)))
    colnames(coverage_flags_effects) <- c("PTE_0.5", "NIE_0.5", "NDE_0.5", "PTE_1", "NIE_1", "NDE_1", "PTE_2", "NIE_2", "NDE_2")
    coverage_flags_effects_marginal <- data.frame(matrix(NA, nrow = n_reps, ncol = 3*length(time_points)))
    colnames(coverage_flags_effects_marginal) <- c("PTE(M)_0.5", "NIE(M)_0.5", "NDE(M)_0.5", "PTE(M)_1", "NIE(M)_1", "NDE(M)_1", "PTE(M)_2", "NIE(M)_2", "NDE(M)_2")
    
    
    # Additional storage for PTE, NIE, and NDE
    effects_df <- data.frame(matrix(NA, nrow = n_reps, ncol = length(time_points) * 3))
    colnames(effects_df) <- c(paste0("PTE_", time_points),
                              paste0("NIE_", time_points),
                              paste0("NDE_", time_points))
    
    effects_sd <- effects_df
    
    # Additional storage for PTE, NIE, and NDE
    effects_df_marginal <- data.frame(matrix(NA, nrow = n_reps, ncol = length(time_points) * 3))
    colnames(effects_df_marginal) <- c(paste0("PTE(M)_", time_points),
                                       paste0("NIE(M)_", time_points),
                                       paste0("NDE(M)_", time_points))
    
    effects_sd_marginal <- effects_df_marginal
    
    # Extract true values for the current scenario
    true_values <- scenarios[[scenario_name]]
    true_values$sigma <- 0.07  # Adding sigma if not already there
    true_effects <- combined_true_typical_effects %>%
      filter(Scenario_SampleSize == paste0(scenario_name, "_", size) & time %in% time_points)
    true_effects_empirical <- combined_true_effects %>%
      filter(Scenario_SampleSize == paste0(scenario_name, "_", size) & time %in% time_points)
    
    # Loop through each repetition
    for (rep in 1:n_reps) {
      # Get the model results for the current repetition, scenario, and sample size
      model_result <- all_results[[paste0("rep_", rep)]][[paste0(scenario_name, "_", size)]]
      
      # Extract summary for parameters
      summary_stats <- model_result[["summary"]][c("beta_os_cov[1]", "gamma", "p", "lambda",
                                                   "beta_trt", "mu_kg[1]",
                                                   "mu_ks[1]", "mu_bsld[1]", "sigma",
                                                   "omega_kg", "omega_ks", "omega_bsld"), ]
      # Extract Rhat values
      rhat_values <- summary_stats$rhat
      
      # Set estimates to NA if Rhat exceeds the threshold
      threshold <- 1.05
      medians_df[rep, ] <- ifelse(rhat_values > threshold, NA, summary_stats$median)
      SDs_df[rep, ] <- ifelse(rhat_values > threshold, NA, summary_stats$sd)
      bias_df[rep, ] <- ifelse(rhat_values > threshold, NA, 
                               summary_stats$median - unlist(true_values[c("beta_os", "gamma", "kappa", "lambda",
                                                                           "beta", "mu_kg", "mu_ks", "mu_M0", "sigma", 
                                                                           "omega_kg", "omega_ks", "omega_M0")]))
      
      # Coverage and Rhat flags
      coverage_flags[rep, ] <- ifelse(rhat_values > threshold, NA,
                                      (unlist(true_values[c("beta_os", "gamma", "kappa", "lambda",
                                                            "beta", "mu_kg", "mu_ks", "mu_M0", "sigma", 
                                                            "omega_kg", "omega_ks", "omega_M0")]) >= summary_stats$q5) &
                                        (unlist(true_values[c("beta_os", "gamma", "kappa", "lambda",
                                                              "beta", "mu_kg", "mu_ks", "mu_M0", "sigma", 
                                                              "omega_kg", "omega_ks", "omega_M0")]) <= summary_stats$q95))
      
      rhat_df[rep, ] <- ifelse(rhat_values > threshold, NA, rhat_values)
      
      # Calculate Rhat flag
      rhat_flags[rep, ] <-  1.05 >= summary_stats$rhat 
      
      if (all(rhat_values <= threshold)) {
        # Extract effects for PTE, NIE, NDE
        effects <- model_result[["effects_typical"]]
        effects_at_time <- effects %>%
          filter(Time %in% c(2, 3, 5)) %>%
          pivot_wider(names_from = Metric, values_from = Median:Upper_95) |>
          mutate(Time = c(0.5, 1.0, 2.0))
        effects_df[rep, ] <- c(effects_at_time$Median_typical_PTE, effects_at_time$Median_typical_NIE, effects_at_time$Median_typical_NDE)
        effects_sd[rep, ] <- c(effects_at_time$SD_typical_PTE, effects_at_time$SD_typical_NIE, effects_at_time$SD_typical_NDE)
        coverage_flags_effects[rep, ] <- true_effects[, c("time", "PTE", "NIE", "NDE")] |> pivot_wider(names_from = time, values_from = PTE:NDE) >=
          effects_at_time[, c("Time", "Lower_95_typical_PTE", "Lower_95_typical_NIE", "Lower_95_typical_NDE")] |> pivot_wider(names_from = Time, values_from = Lower_95_typical_PTE:Lower_95_typical_NDE) &
          true_effects[, c("time", "PTE", "NIE", "NDE")] |> pivot_wider(names_from = time, values_from = PTE:NDE) <=
          effects_at_time[, c("Time", "Upper_95_typical_PTE", "Upper_95_typical_NIE", "Upper_95_typical_NDE")] |> pivot_wider(names_from = Time, values_from = Upper_95_typical_PTE:Upper_95_typical_NDE)
        
        
        # Extract Marginal effects for PTE, NIE, NDE
        effects_marginal <- model_result[["rizopoulos_effects_marginal"]]
        effects_marginal_at_time <- effects_marginal %>%
          filter(Time %in% c(2, 3, 5)) %>%
          pivot_wider(names_from = Metric, values_from = Median:Upper_95) |>
          mutate(Time = c(0.5, 1.0, 2.0))
        effects_df_marginal[rep, ] <- c(effects_marginal_at_time$Median_rizopoulos_marginal_PTE, effects_marginal_at_time$Median_rizopoulos_marginal_NIE, effects_marginal_at_time$Median_rizopoulos_marginal_NDE)
        effects_sd_marginal[rep, ] <- c(effects_marginal_at_time$SD_rizopoulos_marginal_PTE, effects_marginal_at_time$SD_rizopoulos_marginal_NIE, effects_marginal_at_time$SD_rizopoulos_marginal_NDE)
        coverage_flags_effects_marginal[rep, ] <- true_effects_empirical[, c("time", "PTE", "NIE", "NDE")] |> pivot_wider(names_from = time, values_from = PTE:NDE) >=
          effects_marginal_at_time[, c("Time", "Lower_95_rizopoulos_marginal_PTE", "Lower_95_rizopoulos_marginal_NIE", "Lower_95_rizopoulos_marginal_NDE")] |> pivot_wider(names_from = Time, values_from = Lower_95_rizopoulos_marginal_PTE:Lower_95_rizopoulos_marginal_NDE) &
          true_effects_empirical[, c("time", "PTE", "NIE", "NDE")] |> pivot_wider(names_from = time, values_from = PTE:NDE) <=
          effects_marginal_at_time[, c("Time", "Upper_95_rizopoulos_marginal_PTE", "Upper_95_rizopoulos_marginal_NIE", "Upper_95_rizopoulos_marginal_NDE")] |> pivot_wider(names_from = Time, values_from = Upper_95_rizopoulos_marginal_PTE:Upper_95_rizopoulos_marginal_NDE)
      } else {
        # If any Rhat exceeds the threshold, set effects to NA
        effects_df[rep, ] <- NA
        effects_sd[rep, ] <- NA
        effects_df_marginal[rep, ] <- NA
        effects_sd_marginal[rep, ] <- NA
        coverage_flags_effects[rep, ] <- NA
        coverage_flags_effects_marginal[rep, ] <- NA
        
      }
    }
    
    # Calculate mean and SD for PTE, NIE, NDE
    mean_effects <- colMeans(effects_df, na.rm = TRUE)
    sd_effects <- apply(effects_df, 2, sd, na.rm = TRUE)
    true_effects_values <- c(true_effects$PTE, true_effects$NIE, true_effects$NDE)
    effects_bias <- sweep(effects_df, 2, true_effects_values, "-")
    mean_bias_effects <- colMeans(effects_bias, na.rm = TRUE)
    sd_bias_effects <- apply(effects_bias, 2, sd, na.rm = TRUE)
    mean_sd_effects <-  apply(effects_sd, 2, median, na.rm = TRUE) 
    coverage_percentage_effects <- colMeans(coverage_flags_effects, na.rm = TRUE) * 100
    
    # Combine PTE, NIE, and NDE into rows
    effects_table <- data.frame(
      Parameter = colnames(effects_df),
      True_Value = true_effects_values,
      Mean_of_Medians = mean_effects,
      SD_of_Medians = sd_effects,
      Mean_of_SDs = mean_sd_effects,  # Not calculated for effects
      Mean_Bias = mean_bias_effects,
      SD_Bias = sd_bias_effects,
      Relative_Bias_Percent = (mean_bias_effects / true_effects_values) * 100,
      Coverage_Percentage = coverage_percentage_effects,  # Coverage not calculated for effects
      Rhat = NA
    )
    
    ## Marginal effects table 
    # Calculate mean and SD for PTE, NIE, NDE
    mean_effects_marginal <- colMeans(effects_df_marginal, na.rm = TRUE)
    sd_effects_marginal <- apply(effects_df_marginal, 2, sd, na.rm = TRUE)
    true_effects_values <- c(true_effects$PTE, true_effects$NIE, true_effects$NDE)
    effects_bias_marginal <- sweep(effects_df_marginal, 2, true_effects_values, "-")
    mean_bias_effects_marginal <- colMeans(effects_bias_marginal, na.rm = TRUE)
    sd_bias_effects_marginal <- apply(effects_bias_marginal, 2, sd, na.rm = TRUE)
    mean_sd_effects_marginal <-  apply(effects_sd_marginal, 2, median, na.rm = TRUE) 
    coverage_percentage_effects_marginal <- colMeans(coverage_flags_effects_marginal, na.rm = TRUE) * 100
    
    # Combine PTE, NIE, and NDE into rows
    effects_table_marginal <- data.frame(
      Parameter = colnames(effects_df_marginal),
      True_Value = true_effects_values,
      Mean_of_Medians = mean_effects_marginal,
      SD_of_Medians = sd_effects_marginal,
      Mean_of_SDs = mean_sd_effects_marginal,  # Not calculated for effects
      Mean_Bias = mean_bias_effects_marginal,
      SD_Bias = sd_bias_effects_marginal,
      Relative_Bias_Percent = (mean_bias_effects_marginal / true_effects_values) * 100,
      Coverage_Percentage = coverage_percentage_effects_marginal,  # Coverage not calculated for effects
      Rhat = NA
    )
    
    # Calculate the mean of medians and other stats for parameters
    mean_medians <- colMeans(medians_df, na.rm = TRUE)
    sd_medians <- apply(medians_df, 2, sd, na.rm = TRUE)
    mean_of_SDs <- colMeans(SDs_df, na.rm = TRUE)
    mean_rhat <- colMeans(rhat_df, na.rm = TRUE)
    mean_bias <- colMeans(bias_df, na.rm = TRUE)
    sd_bias <- apply(bias_df, 2, sd, na.rm = TRUE)
    relative_bias <- (mean_bias / unlist(true_values[c("beta_os", "gamma", "kappa", "lambda",
                                                       "beta", "mu_kg", "mu_ks", "mu_M0", "sigma", 
                                                       "omega_kg", "omega_ks", "omega_M0")])) * 100
    coverage_percentage <- colMeans(coverage_flags, na.rm = TRUE) * 100
    rhat_percentage <- colSums(rhat_flags, na.rm = TRUE)
    
    # Combine results into a table
    comparison_table <- data.frame(
      Parameter = params_to_extract,
      True_Value = unlist(true_values[c("beta_os", "gamma", "kappa", "lambda",
                                        "beta", "mu_kg", "mu_ks", "mu_M0", "sigma", 
                                        "omega_kg", "omega_ks", "omega_M0")]),
      Mean_of_Medians = mean_medians,
      Mean_Bias = mean_bias,
      SD_Bias = sd_bias,
      Relative_Bias_Percent = relative_bias,
      Coverage_Percentage = coverage_percentage,
      SD_of_Medians = sd_medians,
      Mean_of_SDs = mean_of_SDs,
      Rhat_perc = paste0(round(rhat_percentage / n_reps, 2)*100),
      Rhat = mean_rhat
    )
    
    # Combine parameter and effects tables
    combined_table <- bind_rows(comparison_table, effects_table, effects_table_marginal) |> 
      round_values()
    
    # Add to summary_results list
    summary_results[[paste0(scenario_name, "_", size)]] <- combined_table
  }
}



# Display each table for each scenario and sample size
for (name in names(summary_results)) {
  cat("\nSummary for", name, ":\n")
  print(summary_results[[name]])
}

