
# Load Libraries
require(dplyr)
require(truncnorm)
require(survival)
require(posterior)
require(ggplot2)
require(ggpubr)
require(conflicted)
# we recommend running this in a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
require(cmdstanr)



# Set conflict preferences
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflicts_prefer(posterior::sd)
# ---- Helper Functions ----


# Load custom functions for simulation
source("Functions.R")

stan_model <- cmdstan_model(stan_file = "gq_tgi_beta_ks_os_trt_surv_weibull_link_ks.stan")

sink("package_versions.txt")
sessionInfo()
sink()

# Define parameters for each scenario
scenarios <- list(
  scenario_1 = list(mu_M0 = 60, mu_kg = 0.35, mu_ks = 0.53, beta = 0.5, lambda = 1.2, kappa = 1.5, gamma = -0.2, beta_os = -0.5, omega_ks = 1.3, omega_kg = 1, omega_M0 = .6),
  scenario_2 = list(mu_M0 = 60, mu_kg = 0.35, mu_ks = 0.53, beta = .7, lambda = 1.2, kappa = 1.5, gamma = -1, beta_os = -0.05, omega_ks = 1.3, omega_kg = 1, omega_M0 = .6),
  scenario_3 = list(mu_M0 = 60, mu_kg = 0.35, mu_ks = 0.53, beta = 0.1, lambda = 1.2, kappa = 1.5, gamma = -0.05, beta_os = -0.7, omega_ks = 1.3, omega_kg = 1, omega_M0 = .6)
)

# List of sample sizes for small and large datasets
sample_sizes <- list(small = 200, large = 800)


# Set up parameters
n_reps <- 100
all_results <- list()

# Define parameters to extract
params_to_extract <- c("beta_os_cov", "gamma", "p", "lambda", "beta_trt",
                       "mu_kg[1]", "mu_ks[1]", "mu_bsld", 
                       "omega_kg", "omega_ks", "omega_bsld", "sigma")

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


