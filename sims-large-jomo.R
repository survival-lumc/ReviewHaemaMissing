# Use this to check convergence of the latent categoricals
library(MASS)
library(jomo)
source("jomo-helpers.R")

set.seed(2768)

p <- 6L # number of predictors
covmat <- diag(x = 1, nrow = p + 1L)
n <- 400
betas_true <- rnorm(p + 1L)
dat <- data.frame(mvrnorm(n = 200, mu = betas_true, Sigma = covmat))
colnames(dat) <- c("Y", paste0("X", seq_len(p)))
form <- reformulate(termlabels = paste0("X", seq_len(p)), response = "Y")
form

# Let's make half categorical
p_cat <- p#ceiling(p / 2L)
names_p_cat <- sample(paste0("X", seq_len(p)), replace = FALSE, size = p_cat)


dat[names_p_cat] <- lapply(dat[names_p_cat], function(col) {

  n_cats <- sample(c(2, 3, 4), size = 1, replace = FALSE)
  ggplot2::cut_number(col, n = n_cats)
})

imps_burn <- jomo.lm.MCMCchain(
  formula = form,
  data = dat,
  nburn = 10000,
  out.iter = 100
)

# Test

#profvis::profvis({
#  obj <- assess_jomo_chain(MCMCchain_object = imps_burn, formula_subst = form, dat = dat)
#})

obj <- assess_jomo_chain(MCMCchain_object = imps_burn, formula_subst = form, dat = dat)

# Traces
bayesplot::mcmc_trace(x = obj$betas_subst)
bayesplot::mcmc_trace(x = obj$betas_norm)
bayesplot::mcmc_trace(x = obj$covar_norm)
bayesplot::mcmc_trace(
  x = obj$covar_norm,
  regex_pars = "mdsclass"
)


# ACFs
#lapply(obj, bayesplot::mcmc_acf_bar)
bayesplot::mcmc_acf_bar(x = obj$betas_subst, lags = 500)
bayesplot::mcmc_acf_bar(x = obj$betas_norm, lags = 500)
bayesplot::mcmc_acf_bar(x = obj$covar_norm, lags = 1000, regex_pars = "karnofsk")

# Other
bayesplot::mcmc_dens(x = obj$betas_subst)
bayesplot::mcmc_intervals(x = obj$betas_subst)
bayesplot::mcmc_intervals(x = obj$betas_norm)
bayesplot::mcmc_intervals_data(x = obj$betas_norm)
bayesplot::mcmc_pairs(x = obj$betas_norm)
