prop_slow_models_4PL
models_list <- split_model(
  model = prop_slow_models_4PL,
  split_indices_params = rep(1:J, each = 3),
  split_indices_mu = rep(1:J, each = 2 * (K + 1))
)


  models_list[[1]]$treatment_params_null
