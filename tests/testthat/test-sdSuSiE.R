# =============================================================================
# sdSuSiE() - BASIC FUNCTIONALITY
# =============================================================================
test_that("sdSuSiE returns valid R6 object", {
  #load the simulated data for region 23 in the baseline setting, rs8020953 is the causal variant
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  fit = sdSuSiE(LDmat_list, summary_list)

  expect_true("PIP" %in% names(fit))
  expect_true("PIP_config" %in% names(fit))
  expect_true("alpha" %in% names(fit))
  expect_true("KL" %in% names(fit))
  expect_true("sigma2" %in% names(fit))
  expect_true("V" %in% names(fit))
  expect_true("ELBO" %in% names(fit))
  expect_true("Credible_Set" %in% names(fit))
  expect_true("posterior_mean" %in% names(fit))
  expect_true("posterior_variance" %in% names(fit))
})

test_that("sdSuSiE has correct dimensions", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  fit = sdSuSiE(LDmat_list, summary_list)

  expect_length(fit$PIP, dim(simulated_data$GWAS_M)[1])
  expect_equal(dim(fit$PIP_config), c(dim(simulated_data$GWAS_M)[1], 3))
  expect_length(fit$alpha, 10)
  for(i in 1:10){
  expect_equal(dim(fit$alpha[[i]]), c(dim(simulated_data$GWAS_M)[1], 3))
  }
  expect_length(fit$V, 10)
  for(i in 1:10){
  expect_equal(dim(fit$V[[i]]), c(2, 2))
  }
  expect_equal(dim(fit$posterior_mean), c(dim(simulated_data$GWAS_M)[1], 2))
  expect_equal(dim(fit$posterior_variance), c(dim(simulated_data$GWAS_M)[1], 3))  
})

test_that("sdSuSiE maintains valid probability distributions", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  fit = sdSuSiE(LDmat_list, summary_list)

  # Alpha sum to 1
  for(i in 1:10){
  expect_equal(sum(fit$alpha[[i]]), 1, tolerance = 1e-10)
  }
  # Alpha values are valid probabilities
  for(i in 1:10){
  expect_true(all(fit$alpha[[i]] >= 0 & fit$alpha[[i]] <= 1))
  }
  # PIPs are valid probabilities
  expect_true(all(fit$PIP >= 0 & fit$PIP <= 1))
  expect_true(all(fit$PIP_config >= 0 & fit$PIP_config <= 1))
})

test_that("sdSuSiE ELBO is monotonically increasing", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  fit = sdSuSiE(LDmat_list, summary_list)

  elbo_diff <- na.omit(diff(fit$ELBO))
  expect_true(all(elbo_diff > -1e-6))
})

# =============================================================================
# sdSuSiE() - PARAMETER HANDLING
# =============================================================================

test_that("sdSuSiE respects L parameter", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  fit_L3 = sdSuSiE(LDmat_list, summary_list, L = 3)
  fit_L7 = sdSuSiE(LDmat_list, summary_list, L = 7)

  expect_length(fit_L3$alpha, 3)
  expect_length(fit_L7$alpha, 7)
})

test_that("sdSuSiE handles prior_weights parameter", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  # Uniform weights
  fit_uniform <- sdSuSiE(LDmat_list, summary_list)

  # Custom weights (favor first 10 variables)
  custom_weights <- c(rep(10, 10), rep(1, (dim(simulated_data$GWAS_M)[1]-10)))
  fit_custom <- sdSuSiE(LDmat_list, summary_list, prior_weights = custom_weights)

  expect_true(inherits(fit_uniform, "R6"))
  expect_true(inherits(fit_custom, "R6"))
})

test_that("sdSuSiE handles config_weights parameter", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  # Uniform weights
  fit_uniform <- sdSuSiE(LDmat_list, summary_list)

  # Custom weights (favor the sex-dimorphic effect)
  custom_weights <- c(1/10,7/10,3/10)
  fit_custom <- sdSuSiE(LDmat_list, summary_list, config_weights = custom_weights)

  expect_true(inherits(fit_uniform, "R6"))
  expect_true(inherits(fit_custom, "R6"))
})

test_that("sdSuSiE handles estimate_residual_variance parameter", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  # Uniform weights
  fit_uniform <- sdSuSiE(LDmat_list, summary_list)

  # Estimated the residual variance
  fit_estimated <- sdSuSiE(LDmat_list, summary_list, estimate_residual_variance  = TRUE)

  expect_true(inherits(fit_uniform, "R6"))
  expect_true(inherits(fit_estimated, "R6"))
})

# =============================================================================
# sdSuSiE() - SIGNAL RECOVERY
# =============================================================================

test_that("sdSuSiE identifies true causal variables", {
  data(simulated_data)

  summary_list = list()
  summary_list[[1]] = simulated_data$GWAS_M
  summary_list[[2]] = simulated_data$GWAS_F
  names(summary_list) = c("Male","Female")

  LDmat_list = list()
  LDmat_list[[1]] = simulated_data$LD_M
  LDmat_list[[2]] = simulated_data$LD_F
  names(LDmat_list) = c("Male","Female")

  fit = sdSuSiE(LDmat_list, summary_list)

  # Top PIPs should include true causal variables
  top_vars <- GWAS_M[order(fit$PIP, decreasing = TRUE),]$SNP[1:10]
  overlap <- length(intersect(top_vars, "rs8020953"))

  expect_true(overlap >= 1)
})