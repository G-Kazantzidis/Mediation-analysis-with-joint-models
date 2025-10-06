functions { 
  //- ifelse ----
  // Vectorized ifelse() similar as in R.
  
  row_vector ifelse(array[] int condition, row_vector yes, row_vector no) {
    row_vector[num_elements(yes)] result;
    for (i in 1:num_elements(yes)) {
      result[i] = condition[i] ? yes[i] : no[i];
    }
    return result;
  }
  
  //- is_negative ----
  // Vectorized negative predicate. Basically returns ifelse(x < 0, 1, 0).
  
  array[] int is_negative(row_vector x) {
    array[num_elements(x)] int result;
    for (i in 1:num_elements(x)) {
      result[i] = x[i] < 0.0;
    }
    return result;
  }
  
  
  //- row_means ----
  // Means across the rows of a matrix
  row_vector row_means(matrix x) {
    row_vector[cols(x)] result = rep_row_vector(1.0 / rows(x), rows(x)) * x;
    return result;
  }
  
  //- is_positive ----
  // Vectorized positive predicate. Basically returns ifelse(x > 0, 1, 0).

  array[] int is_positive(row_vector x) {
    return is_negative(- x);
  }
  
  //- which ----
  // which function. Returns which of the 0/1 elements are not 0.

  array[] int which(array[] int x) {
    int len = sum(x);
    array[len] int ret;
    int pointer = 0;
    for (i in 1:num_elements(x)) {
      if (x[i] != 0) {
        if (x[i] != 1) {
          reject("integer array passed to `which` function must only contain 0 or 1");
        }
        pointer = pointer + 1;
        ret[pointer] = i;
      }
    }
    return ret;
  }
  
  // SLD model 
  row_vector sld(row_vector time, row_vector psi_bsld, row_vector psi_ks, 
  row_vector psi_kg, row_vector trt, real beta_trt) {
    row_vector[num_elements(time)] result = fmin(
      3000.0,
      ifelse(is_negative(time), 
      psi_bsld .* ( exp(psi_kg .* time)), 
      psi_bsld .* ( exp(psi_kg .* time) + exp(-(psi_ks .* (exp(beta_trt .* trt))) .* time) - 1)));
    return result;
  }
  
  // SLD model per visit
  real sld_per_visit(real time, 
  real psi_bsld, real psi_ks, 
  real psi_kg, real trt, real beta_trt) {
    real result = psi_bsld .* ( exp(psi_kg .* time) + exp(-(psi_ks .* (exp(beta_trt .* trt))) .* time) - 1);
    return result;
  }
  
  // SLD model per visit
  real sld_neg_visit(real time, 
  real psi_bsld, real psi_kg) {
    real result = psi_bsld .* ( exp(psi_kg .* time));
    return result;
  }
  
  //- neg_log_sqrt_2_pi ----
  // Constant used in below.
  real neg_log_sqrt_2_pi() {
    return -0.9189385332046727;
  }
  
  //- vect_normal_log_dens ----
  // Vectorized version of normal_lpdf, i.e. returns log normal density values.
  row_vector vect_normal_log_dens(row_vector y, row_vector mu, row_vector sigma) {
    row_vector[num_elements(y)] y_stand = (y - mu) ./ sigma;
    row_vector[num_elements(y)] main_result = - (y_stand .* y_stand) / 2;
    return main_result + neg_log_sqrt_2_pi() - log(sigma);
  }
  
    //- vect_normal_log_cum ----
  // Vectorized version of normal_lcdf, i.e. return log normal CDF values.
  row_vector vect_normal_log_cum(real quantile, row_vector mu, row_vector sigma) {
    row_vector[num_elements(mu)] quant_stand = (quantile - mu) ./ sigma;
    row_vector[num_elements(mu)] cdf_vals = Phi(quant_stand);
    return log(cdf_vals);
  }
  
  
  
  //- curr_sld ----
  // current SLD value 

  matrix curr_sld(matrix time, row_vector psi_bsld, row_vector psi_ks,
  row_vector psi_kg, row_vector trt_pts, row_vector psi_trt) {
    // Here we assume that psi's are replicated along the rows of the time matrix.
    matrix[rows(time), cols(psi_bsld)] psi_bsld_matrix = rep_matrix(psi_bsld, rows(time));
    matrix[rows(time), cols(psi_ks)] psi_ks_matrix = rep_matrix(psi_ks, rows(time));
    matrix[rows(time), cols(psi_kg)] psi_kg_matrix = rep_matrix(psi_kg, rows(time));
    matrix[rows(time), cols(psi_kg)] psi_trt_matrix = rep_matrix(psi_trt, rows(time));
    matrix[rows(time), cols(psi_kg)] trt_pts_mat = rep_matrix(trt_pts, rows(time));

    // We also assume that all the time values are positive. 
    matrix[rows(time), cols(time)] result = fmin(
      3000.0,
      psi_bsld_matrix .* (exp(psi_kg_matrix .* time) + exp(-(psi_ks_matrix .* (exp(trt_pts_mat .* psi_trt_matrix))).* time ) - 1));
      
    return result;
  }

  
 //- log_h0 ----
  // Baseline log hazard
  
  matrix log_h0(matrix time, real lambda, real p) {
  matrix[rows(time), cols(time)] result;
  
  // Weibull structure
  result = log(p) - log(lambda) + (p - 1) * log(time ./ lambda);
  
  return result;
}

  //- log_hazard ----
  // Log Hazard function

matrix log_hazard(matrix time,
                  real lambda, real p, real gamma,
                  row_vector eta_tilde_ks, real beta_trt,
                  vector beta_os_cov, matrix os_cov_design, matrix os_cov_design_tgi) {
  matrix[rows(time), cols(time)] log_baseline = log_h0(time, lambda, p);

  // Convert eta_tilde_ks to vector for element-wise multiplication
  vector[rows(os_cov_design_tgi)] eta_tilde_ks_vec = to_vector(log(eta_tilde_ks));
  
    // Calculate the treatment effect
    vector[rows(os_cov_design_tgi)] treatment_effect = os_cov_design_tgi[, 1] * beta_trt;
    // print("Dimensions of treatment_effect: ", dims(treatment_effect));
    // print("First 20 elements of treatment_effect: ", treatment_effect[1:10]);
    
    
    // Element-wise multiply eta_tilde_kg with treatment_effect
    vector[rows(os_cov_design_tgi)] eta_tilde_ks_effect = eta_tilde_ks_vec + treatment_effect ; // Add here the term + treatment_effect 
    // print("Dimensions of eta_tilde_kg_effect: ", dims(eta_tilde_kg_effect));
    // print("First 20 elements of eta_tilde_kg_effect: ", eta_tilde_kg_effect[1:10]);


  // Replicate the result into a matrix
  matrix[rows(time), cols(eta_tilde_ks)] eta_tilde_ks_matrix = rep_matrix(eta_tilde_ks_effect', rows(time));

  row_vector[rows(os_cov_design)] cov_contribution = (os_cov_design * beta_os_cov)';
  matrix[rows(time), cols(cov_contribution)] cov_contribution_matrix = rep_matrix(cov_contribution, rows(time));

  matrix[rows(time), cols(time)] result = log_baseline + gamma * eta_tilde_ks_matrix + cov_contribution_matrix;
  return result;
}

  
  //- log_survival ----
  // Survival function
  // now we take care of time = 0 values - for these we just directly know that the result is 0.

  row_vector log_survival(row_vector time, real lambda,
                        real p, real gamma,
                        row_vector eta_tilde_ks, real beta_trt,
                        data vector nodes, data row_vector weights, vector beta_os_cov, 
                        data matrix os_cov_design, data matrix os_cov_design_tgi) {
  array[cols(time)] int time_positive = is_positive(time);
  int n_positive = sum(time_positive);
  array[n_positive] int time_positive_index = which(time_positive);
  row_vector[cols(time)] result = rep_row_vector(0.0, cols(time));

  matrix[rows(nodes), n_positive] nodes_time = (nodes + 1) * (time[time_positive_index] / 2);
  matrix[rows(nodes), n_positive] nodes_time_hazard = fmin(8000.0, exp(log_hazard(nodes_time, lambda,
                                                                                   p, gamma,
                                                                                   eta_tilde_ks[time_positive_index], beta_trt,
                                                                                   beta_os_cov, os_cov_design[time_positive_index], 
                                                                                   os_cov_design_tgi[time_positive_index])));
  
  result[time_positive_index] = - (weights * nodes_time_hazard) .* time[time_positive_index] / 2;
  return result;
}
  
    
  // Cumulative hazard function for survival data
  row_vector cumulative_hazard(row_vector time, real lambda, real p, real gamma,
  row_vector psi_ks, real beta_trt, data vector nodes, data row_vector weights,
  vector beta_os_cov, data matrix os_cov_design, matrix os_cov_design_tgi) {
    
    array[cols(time)] int time_positive = is_positive(time);
    int n_positive = sum(time_positive);
    array[n_positive] int time_positive_index = which(time_positive);
    row_vector[cols(time)] result = rep_row_vector(0.0, cols(time));
    
    if (n_positive > 0) {
      matrix[rows(nodes), n_positive] nodes_time = (nodes + 1) * (time[time_positive_index] / 2);
      matrix[rows(nodes), n_positive] nodes_time_hazard = fmin(8000.0, exp(log_hazard(nodes_time, lambda, p, gamma,
      psi_ks[time_positive_index], beta_trt, 
      beta_os_cov, os_cov_design[time_positive_index], os_cov_design_tgi[time_positive_index]
      )));
      result[time_positive_index] = (weights * nodes_time_hazard) .* time[time_positive_index] / 2;
    }
    return result;
  }
  
  
}

data {
  // Longitudinal data 
  int<lower=1> Nind; // Number of individuals.
  int<lower=1> Nind_cens; // Number of individuals with censored sld values.
  int<lower=1> Nta_total; // Total number of tumor assessments.
  int<lower=1> Nind_dead; // Number of dead individuals (observed survival time).
  int<lower=1> Nta_cens_y; // Number of censored tumor assessments (below threshold).
  int<lower=1> Nta_obs_y; // Number of censored tumor assessments (below threshold).

  int<lower=1> n_arms; // Total number of arms.
  array[Nind] int<lower=1, upper=n_arms> arm_index; // Index of treatment arm for all individuals.
  array[Nind] int<lower=1, upper=2> arm_index_real;

  
  array[Nta_total] int ind_index; // Index of individuals for each tumor assessment.
  array[Nind] int os_ind_index; // Index of individuals for each tumor assessment.
  array[Nta_obs_y] int obs_y_index; // Index of tumor assessments .
  array[Nta_cens_y] int cens_y_index; // Index of censored tumor assessments.
  array[Nind_dead] int dead_ind_index; // Index of dead individuals (observed survival time).
  
  row_vector[Nta_total] Yobs; // Array of individual responses.
  row_vector[Nta_total] Tobs; // Individual timepoints.
  row_vector[Nind] Tobs_neg; // Individual timepoints.
  row_vector[Nta_total] trt; // trt per patient  
  
  
   real Ythreshold; // Censoring threshold.
  // Matrix of individuals x observed tumor assessments (sparse matrix of 0s and 1s),
  // so the dimension is Nind x Nta_obs_y.
  int<lower=1> n_w_mat_inds_obs_y;
  vector[n_w_mat_inds_obs_y] w_mat_inds_obs_y;
  int<lower=1> n_v_mat_inds_obs_y;
  array[n_v_mat_inds_obs_y] int v_mat_inds_obs_y;
  int<lower=1> n_u_mat_inds_obs_y;
  array[n_u_mat_inds_obs_y] int u_mat_inds_obs_y;
  
  // Matrix of individuals x censored tumor assessments (sparse matrix of 0s and 1s).
  // so the dimension is Nind x Nta_cens_y.
  int<lower=1> n_w_mat_inds_cens_y;
  vector[n_w_mat_inds_cens_y] w_mat_inds_cens_y;
  int<lower=1> n_v_mat_inds_cens_y;
  array[n_v_mat_inds_cens_y] int v_mat_inds_cens_y;
  int<lower=1> n_u_mat_inds_cens_y;
  array[n_u_mat_inds_cens_y] int u_mat_inds_cens_y;

  // Matrix of individuals x observed tumor assessments (sparse matrix of 0s and 1s),
  // so the dimension is Nind x Nta_obs_y.
  int<lower=1> n_w_mat_inds_obs_y_all;
  vector[n_w_mat_inds_obs_y_all] w_mat_inds_obs_y_all;
  int<lower=1> n_v_mat_inds_obs_y_all;
  array[n_v_mat_inds_obs_y_all] int v_mat_inds_obs_y_all;
  int<lower=1> n_u_mat_inds_obs_y_all;
  array[n_u_mat_inds_obs_y_all] int u_mat_inds_obs_y_all;
  
  // prediction data  
  row_vector[Nind] final_sld_time; // final sld obs time 
  row_vector[Nind] final_sld;  // final sld obs 
  
  row_vector[Nind] trt_gen_q;  // treatment per patient
    
  // Times for SLD predictions.
  int<lower=1> n_sld_exp_times;
  int<lower=1> n_sld_pred_times;
  
  row_vector<lower=0>[n_sld_exp_times] sld_exp_times;
  row_vector<lower=0>[n_sld_exp_times] sld_pred_times;
  int<lower=1> total_obs_sld_times;
  
  // Survival data ----
  row_vector[Nind] Times;
  row_vector[Nind] trt_pts;
  int<lower=1> p_os_cov_design;
  matrix[Nind, p_os_cov_design] os_cov_design;
  matrix[Nind, p_os_cov_design] s1_design;
  matrix[Nind, p_os_cov_design] s0_design;
  
  // Times for OS predictions.
  int<lower=1> n_os_pred_times;
  row_vector<lower=0>[n_os_pred_times] os_pred_times;
  
  // Integration parameters ----
  int<lower=1> n_nodes;
  vector[n_nodes] nodes;
  row_vector<lower=0, upper=1>[n_nodes] weights;
}


parameters {
  
     // Hyper parameters
  real<lower=0, upper=100> mean_mu_ks;
  real<lower=0, upper=100> mean_mu_kg;
  real<lower=0, upper=100> mean_mu_bsld;
  
  real<lower=0, upper=10> sd_mu_ks;
  real<lower=0, upper=10> sd_mu_kg;
  real<lower=0, upper=10> sd_mu_bsld;
  
  // Standard deviation for RUV.
  real<lower=0.00001, upper=100> sigma;

  
  // Population parameters.
  row_vector<lower=0>[n_arms] mu_bsld;
  row_vector<lower=0, upper=100>[n_arms] mu_ks;
  row_vector<lower=0, upper=100>[n_arms] mu_kg;
  
  real<lower=0> omega_bsld;
  real<lower=0> omega_ks;
  real<lower=0> omega_kg;
  
  // Random effects.
  row_vector[Nind] eta_tilde_bsld;
  row_vector[Nind] eta_tilde_ks;
  row_vector[Nind] eta_tilde_kg;
  
  real<lower=-10> beta_trt;
  
  // Survival parameters.
  real p;
  real<lower=0> lambda; // For the log-logistic baseline hazard..
  real gamma; // Link parameter for TTG.
  vector[p_os_cov_design] beta_os_cov; // Covariate coefficients.
}


transformed parameters {
  // Non-centered reparametrization for hierarchical models.
  row_vector[Nind] psi_bsld = exp(log(mu_bsld[arm_index]) + eta_tilde_bsld * omega_bsld); 
  row_vector[Nind] psi_ks = exp(log(mu_ks[arm_index]) + eta_tilde_ks * omega_ks);
  row_vector[Nind] psi_kg = exp(log(mu_kg[arm_index]) + eta_tilde_kg * omega_kg);
  row_vector[Nind] psi_trt = rep_row_vector(beta_trt, Nind);  
  
  // Log-likelihood values for using the loo package.
  row_vector[Nind] log_lik;
  row_vector[Nta_total] Ypred;
  row_vector[Nta_total] Ypred_log_lik;
  row_vector[Nind] log_surv_vals;
    
  log_lik = rep_row_vector(0.0, Nind);
  
  
  Ypred = sld(Tobs,
  psi_bsld[ind_index], psi_ks[ind_index], psi_kg[ind_index], trt, beta_trt);
  
  Ypred_log_lik[obs_y_index] = vect_normal_log_dens(Yobs[obs_y_index], Ypred[obs_y_index], Ypred[obs_y_index] * sigma);
  Ypred_log_lik[cens_y_index] = vect_normal_log_cum(Ythreshold, Ypred[cens_y_index], Ypred[cens_y_index] * sigma);
  
  log_lik += csr_matrix_times_vector(Nind, Nta_total, w_mat_inds_obs_y_all, v_mat_inds_obs_y_all, u_mat_inds_obs_y_all,
  Ypred_log_lik')';
  
  // We always add the log-survival to the log-likelihood.
  log_surv_vals = log_survival(Times[os_ind_index], lambda, p, gamma,
  psi_ks[os_ind_index], beta_trt, 
  nodes, weights, beta_os_cov, os_cov_design[os_ind_index], os_cov_design[os_ind_index]);
  log_lik += log_surv_vals;
  
   // In case of death we add the log-hazard on top.
  log_lik[dead_ind_index] += to_row_vector(log_hazard(to_matrix(Times[dead_ind_index]), lambda, p, gamma,
  psi_ks[dead_ind_index], beta_trt,
  beta_os_cov, os_cov_design[dead_ind_index], os_cov_design[dead_ind_index]));
  
}

model {

  mean_mu_bsld ~ lognormal(1, 1);
  
  mean_mu_ks ~ normal(1,1); 
  mean_mu_kg ~ normal(0, .5); 
  
  sd_mu_ks ~ lognormal(0,0.5);
  sd_mu_kg ~ lognormal(0,0.5);
  sd_mu_bsld ~ lognormal(1,0.5);
  
  mu_ks[n_arms] ~ lognormal(mean_mu_ks, sd_mu_ks);
  mu_kg[n_arms] ~ lognormal(mean_mu_kg, sd_mu_kg);
  mu_bsld[n_arms] ~ lognormal(mean_mu_bsld, sd_mu_bsld); 
  
  sigma ~ lognormal(-1.6, 0.8); // Proportional error
    
  omega_bsld ~ lognormal(0,1);
  omega_ks ~ lognormal(0,1);
  omega_kg ~ lognormal(0,1);
  
  eta_tilde_bsld ~ normal(0,1);
  eta_tilde_ks ~ normal(0,1);
  eta_tilde_kg ~ normal(0,1);
  
  p ~ lognormal(0, 0.5);
  lambda ~ lognormal(0, 0.5);
  gamma ~ normal(0, 2);
  beta_os_cov ~ normal(0, 5); 

  
  target += sum(log_lik);
}


generated quantities {
  matrix[Nind, n_sld_exp_times] ind_expected_sld; // Timepoints to estimate SLD 
  matrix[Nind, n_sld_pred_times] ind_predicted_sld; // Timepoints to predict SLD 
  matrix[Nind, n_sld_exp_times] ind_exp_treated; // Timepoints to estimate SLD 
  matrix[Nind, n_sld_pred_times] ind_exp_not_treated; // Timepoints to predict SLD 
  matrix[Nind, n_sld_pred_times] ind_exp_control; // Timepoints to predict SLD 
  array[n_sld_exp_times] int i_rep;  // empty area for estimated parameters 
  row_vector[Nind] treated = rep_row_vector(1, Nind);
  row_vector[Nind] not_treated = rep_row_vector(0, Nind);
  row_vector[Nind] zero_ks = rep_row_vector(0, Nind);
  
  row_vector[n_os_pred_times] pop_PTE;
  
  matrix[Nind, n_os_pred_times] s00;
  matrix[Nind, n_os_pred_times] s01;
  matrix[Nind, n_os_pred_times] s000;
  matrix[Nind, n_os_pred_times] s010;
  matrix[Nind, n_os_pred_times] s10;
  matrix[Nind, n_os_pred_times] s11;
  matrix[Nind, n_os_pred_times] NIE;
  matrix[Nind, n_os_pred_times] NDE;
  matrix[Nind, n_os_pred_times] PTE;
  matrix[Nind, n_os_pred_times] NIE_full;
  matrix[Nind, n_os_pred_times] NDE_full;
  matrix[Nind, n_os_pred_times] PTE_full;
  
  row_vector[total_obs_sld_times] ind_expected_sld_exact_times; // Timepoints to estimate SLD at the observed times
  row_vector[Nind] ind_expected_sld_exact_negative_times; // Timepoints to estimate SLD at the observed times 
  
  // Process each patient.
  for (i in 1:Nind) {
    i_rep = rep_array(i, n_sld_exp_times);
    
    // Expected SLD values.
    ind_expected_sld[i] = sld(sld_exp_times,
    psi_bsld[i_rep], psi_ks[i_rep], psi_kg[i_rep], trt_gen_q[i_rep], beta_trt);
    
    // Predicted SLD values.
    ind_predicted_sld[i] = sld(final_sld_time[i_rep] + sld_pred_times,
    final_sld[i_rep], psi_ks[i_rep], psi_kg[i_rep], trt_gen_q[i_rep], beta_trt);
    
    s01[i] = exp(log_survival(os_pred_times,
    lambda, p, gamma,
    psi_ks[i_rep], beta_trt,
    nodes, weights, beta_os_cov, s1_design[i_rep], s0_design[i_rep]));
    
    // Survival not Treated, SLD not treated
    s00[i] = exp(log_survival(os_pred_times,
    lambda, p, gamma,
    psi_ks[i_rep], beta_trt,
    nodes, weights, beta_os_cov, s0_design[i_rep], s0_design[i_rep]));
    
    // Survival Treated, SLD treated
    s11[i] = exp(log_survival(os_pred_times,
    lambda, p, gamma,
    psi_ks[i_rep], beta_trt,
    nodes, weights, beta_os_cov, s1_design[i_rep], s1_design[i_rep]));
    
    // Survival not Treated, SLD treated
    s10[i] = exp(log_survival(os_pred_times,
    lambda, p, gamma,
    psi_ks[i_rep], beta_trt,
    nodes, weights, beta_os_cov, s0_design[i_rep], s1_design[i_rep]));
    
    // NIE
    NIE[i] = s11[i] - s01[i];
    
    // NDE
    NDE[i] = s01[i] - s00[i];
    
    PTE[i] = (s11[i] - s01[i]) ./ (s11[i] - s00[i]);
    
    
  }
  
  
   // Arrays to store posterior samples for population-level metrics
  vector[n_sld_exp_times] NIE_pop_samples;
  vector[n_sld_exp_times] NDE_pop_samples;
  vector[n_sld_exp_times] PTE_pop_samples;

  // // Calculate population-level metrics
  // for (t in 1:n_sld_exp_times) {
  //   // Mean survival probabilities across individuals at time t
  //   real mean_NIE = mean(to_vector(NIE[:, t]));
  //   real mean_NDE = mean(to_vector(NDE[:, t]));
  //   real mean_PTE = mean(to_vector(PTE[:, t]));
  // 
  //   // NIE, NDE, and PTE calculations
  //   NIE_pop_samples[t] = mean_NIE;
  //   NDE_pop_samples[t] = mean_NDE;
  //   PTE_pop_samples[t] = mean_PTE;
  // }
  // Calculate population-level metrics excluding NaN values
for (t in 1:n_sld_exp_times) {
  real sum_NIE = 0.0;
  real sum_NDE = 0.0;
  real sum_PTE = 0.0;
  int count_NIE = 0;
  int count_NDE = 0;
  int count_PTE = 0;

  for (i in 1:Nind) {
    // Check for valid NIE values and accumulate
    if (!is_nan(NIE[i, t])) {
      sum_NIE += NIE[i, t];
      count_NIE += 1;
    }
    // Check for valid NDE values and accumulate
    if (!is_nan(NDE[i, t])) {
      sum_NDE += NDE[i, t];
      count_NDE += 1;
    }
    // Check for valid PTE values and accumulate
    if (!is_nan(PTE[i, t])) {
      sum_PTE += PTE[i, t];
      count_PTE += 1;
    }
  }

  // Calculate means only if there are valid entries
  NIE_pop_samples[t] = (count_NIE > 0) ? sum_NIE / count_NIE : 0.0;
  NDE_pop_samples[t] = (count_NDE > 0) ? sum_NDE / count_NDE : 0.0;
  PTE_pop_samples[t] = (count_PTE > 0) ? sum_PTE / count_PTE : 0.0;
}

  
    matrix[n_os_pred_times, 4] marginal_survival; // Stores s00, s01, s10, s11
  row_vector[n_os_pred_times] typical_NIE;     // Marginal NIE over time
  row_vector[n_os_pred_times] typical_NDE;     // Marginal NDE over time
  row_vector[n_os_pred_times] typical_PTE;     // Marginal PTE over time

  // Calculate marginal survival probabilities at each time point

// Marginal s00: No treatment, no mediator effect
marginal_survival[:, 1] = to_vector( exp(log_survival(os_pred_times, 
  lambda, p, gamma,
  to_row_vector(rep_array(mu_ks[1], n_os_pred_times)), beta_trt,
  nodes, weights, beta_os_cov, s0_design, s0_design)));

// Marginal s01: No treatment, mediator effect
marginal_survival[:, 2] = to_vector(exp(log_survival(os_pred_times, 
  lambda, p, gamma,
  to_row_vector(rep_array(mu_ks[1], n_os_pred_times)), beta_trt,
  nodes, weights, beta_os_cov, s1_design, s0_design)));

// Marginal s10: Treatment, no mediator effect
marginal_survival[:, 3] = to_vector(exp(log_survival(os_pred_times, 
  lambda, p, gamma,
  to_row_vector(rep_array(mu_ks[1], n_os_pred_times)), beta_trt,
  nodes, weights, beta_os_cov, s0_design, s1_design)));

// Marginal s11: Treatment, mediator effect
marginal_survival[:, 4] = to_vector(exp(log_survival(os_pred_times, 
  lambda, p, gamma,
  to_row_vector(rep_array(mu_ks[1], n_os_pred_times)), beta_trt,
  nodes, weights, beta_os_cov, s1_design, s1_design)));




  // Compute marginal NIE, NDE, and PTE over time
  for (t in 1:n_os_pred_times) {
    typical_NIE[t] = marginal_survival[t, 4] - marginal_survival[t, 2];
    typical_NDE[t] = marginal_survival[t, 2] - marginal_survival[t, 1];
    typical_PTE[t] = typical_NIE[t] / (marginal_survival[t, 4] - marginal_survival[t, 1]);
  }
  
  
  
   matrix[Nind, n_os_pred_times] NIE_samples;  // Store NIE for each subject and sample
    matrix[Nind, n_os_pred_times] NDE_samples;  // Store NDE for each subject and sample
    matrix[Nind, n_os_pred_times] PTE_samples;  // Store PTE for each subject and sample
    matrix[100, n_os_pred_times] s00_marginal;
    matrix[100, n_os_pred_times] s01_marginal;
    matrix[100, n_os_pred_times] s10_marginal;
    matrix[100, n_os_pred_times] s11_marginal;
    matrix[100, n_os_pred_times] NIE_marginal;
    matrix[100, n_os_pred_times] NDE_marginal;
    matrix[100, n_os_pred_times] PTE_marginal;
    
    
    row_vector[n_os_pred_times] rizopoulos_marginal_NIE;
    row_vector[n_os_pred_times] rizopoulos_marginal_NDE;
    row_vector[n_os_pred_times] rizopoulos_marginal_PTE;
    
    // Average effects over Monte Carlo samples for each subject
    matrix[Nind, n_os_pred_times] NIE_averaged;
    matrix[Nind, n_os_pred_times] NDE_averaged;
    matrix[Nind, n_os_pred_times] PTE_averaged;
    
    // Monte Carlo integration for each subject
    for (i in 1:Nind) {
      
      
      i_rep = rep_array(i, n_sld_exp_times);
      
      for (s in 1:100) {
        // Sample random effects for subject i
        real sampled_eta_ks = normal_rng(0, omega_ks);

        // Compute individual parameters based on sampled random effects
        real psi_ks_sampled = exp(log(mu_ks[arm_index[i]]) + sampled_eta_ks);

        // Compute survival probabilities for this sample
        s00_marginal[s,] = exp(log_survival(os_pred_times, lambda, p, gamma, 
        to_row_vector(rep_array(psi_ks_sampled, n_os_pred_times)), beta_trt,
        nodes, weights, beta_os_cov, s0_design[i_rep], s0_design[i_rep]));
        
        s01_marginal[s,] = exp(log_survival(os_pred_times, lambda, p, gamma, 
        to_row_vector(rep_array(psi_ks_sampled, n_os_pred_times)), beta_trt,
        nodes, weights, beta_os_cov, s1_design[i_rep], s0_design[i_rep]));
        
        s10_marginal[s,] = exp(log_survival(os_pred_times, lambda, p, gamma, 
        to_row_vector(rep_array(psi_ks_sampled, n_os_pred_times)), beta_trt,
        nodes, weights, beta_os_cov, s0_design[i_rep], s1_design[i_rep]));
        
        s11_marginal[s,] = exp(log_survival(os_pred_times, lambda, p, gamma, 
        to_row_vector(rep_array(psi_ks_sampled, n_os_pred_times)), beta_trt,
        nodes, weights, beta_os_cov, s1_design[i_rep], s1_design[i_rep]));
        
        // Compute NIE, NDE, and PTE for this sample
        NIE_marginal[s, ] = s11_marginal[s, ] - s01_marginal[s, ];
        NDE_marginal[s, ] = s01_marginal[s, ] - s00_marginal[s, ];
        for (t in 1:n_sld_exp_times) {
          PTE_marginal[s, t] = NIE_marginal[s, t] / (s11_marginal[s, t] - s00_marginal[s, t]);
        }
      }

      for (t in 1:n_sld_exp_times) {
        NIE_averaged[i, t] = mean(NIE_marginal[, t]);
        NDE_averaged[i, t] = mean(NDE_marginal[, t]);
        
        // Exclude NaN PTE values during averaging
        real sum_PTE = 0.0;
        int count_PTE = 0;
        
        for (s in 1:100) {
          if (!is_nan(PTE_marginal[s, t])) {
            sum_PTE += PTE_marginal[s, t];
            count_PTE += 1;
          }
        }
        PTE_averaged[i, t] = (count_PTE > 0) ? sum_PTE / count_PTE : 0.0;
      }
    }
    
    // Compute population-level marginal effects
    for (t in 1:n_os_pred_times) {
      rizopoulos_marginal_NIE[t] = mean(NIE_averaged[, t]);
      rizopoulos_marginal_NDE[t] = mean(NDE_averaged[, t]);
      rizopoulos_marginal_PTE[t] = mean(PTE_averaged[, t]);
    }
  // 
  // 
  // // Expected SLD values at exact times
  // for (t in 1:total_obs_sld_times) {
  //   ind_expected_sld_exact_times[t] = sld_per_visit(Tobs[t],
  //   psi_bsld[ind_index[t]], psi_ks[ind_index[t]], psi_kg[ind_index[t]], trt_gen_q[ind_index[t]], beta_trt);
  // }
  // 
  // for (t in 1:Nind) {
  //   ind_expected_sld_exact_negative_times[t] = sld_neg_visit(Tobs_neg[t],
  //   psi_bsld[ind_index[t]], psi_kg[ind_index[t]]);
  // }
  // 
  // matrix[Nind, n_os_pred_times] survival_prob;
  // row_vector[n_os_pred_times] mean_survival_prob;
  // 
  // for (i in 1:Nind) {
  //   for (j in 1:n_os_pred_times) {
  //     survival_prob[i, j] = exp(log_survival(to_row_vector(rep_array(os_pred_times[j], 1)), lambda, p, gamma,
  //     to_row_vector(rep_array(psi_ks[i], 1)), beta_trt,
  //     nodes, weights, beta_os_cov,
  //     to_matrix(os_cov_design[i, ]),
  //     to_matrix(os_cov_design[i, ]))[1]);
  //   }
  // }
  // 
  // for (j in 1:n_os_pred_times) {
  //   mean_survival_prob[j] = mean(survival_prob[, j]);
  // }
  // 
  // matrix[2, n_os_pred_times] mean_survival_prob_arm = rep_matrix(0.0, 2, n_os_pred_times);
  // for (a in 1:2) {
  //   for (j in 1:n_os_pred_times) {
  //     real sum_surv_prob = 0;
  //     int count = 0;
  //     for (i in 1:Nind) {
  //       if (arm_index_real[i] == a) {
  //         sum_surv_prob += survival_prob[i, j];
  //         count += 1;
  //       }
  //     }
  //     mean_survival_prob_arm[a, j] = sum_surv_prob / count;
  //   }
  // }
  // 
  // // New generated quantities for martingale residuals
  // row_vector[Nind] martingale_residuals;
  // 
  // // Process each patient to calculate martingale residuals
  // for (i in 1:Nind) {
  //   // Calculate the cumulative hazard up to the observed time
  //   row_vector[1] obs_time = to_row_vector(rep_array(Times[i], 1));
  //   row_vector[1] cumulative_hazard_val = cumulative_hazard(obs_time, lambda, p, gamma,
  //   to_row_vector(rep_array(psi_ks[i], 1)), 
  //   beta_trt, nodes, weights,
  //   beta_os_cov, to_matrix(os_cov_design[i, ]), to_matrix(os_cov_design[i, ]));
  //   
  //   // Check if the individual is in the dead_ind_index
  //   int event = 0;
  //   for (j in 1:Nind_dead) {
  //     if (dead_ind_index[j] == i) {
  //       event = 1;
  //       break;
  //     }
  //   }
  //   
  //   // Calculate the martingale residual for the individual
  //   martingale_residuals[i] = event - cumulative_hazard_val[1];
  // }
}


