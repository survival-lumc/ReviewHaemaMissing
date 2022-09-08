prepare_param_dfs <- function(MCMCchain_object,
                              formula_subst,
                              dat,
                              intercept = FALSE) {

  # Get labels for subst model
  terms_subst <- dimnames(model.matrix.default(object = formula_subst, data = dat))[[2]]
  if (!intercept) terms_subst <- terms_subst[-1L]

  # Start with betas (y and substantive model)
  # Change this bit with array subsetting (see untitled)
  betas_norm <- as.data.frame.table(MCMCchain_object[["collectbeta"]])
  betas_subst <- as.data.frame.table(MCMCchain_object[["collectbetaY"]])
  betas_norm[["Var1"]] <- betas_subst[["Var1"]] <- NULL
  colnames(betas_norm) <- colnames(betas_subst) <- c("term", "iter", "value")
  betas_subst[["iter"]] <- as.numeric(betas_subst[["iter"]])
  betas_norm[["iter"]] <- as.numeric(betas_norm[["iter"]])
  levels(betas_subst[["term"]]) <- terms_subst

  # Now covariance matrices
  omega_norm <- as.data.frame.table(MCMCchain_object[["collectomega"]])
  colnames(omega_norm) <- c("Var1", "Var2", "iter", "value")
  omega_norm[["iter"]] <- as.numeric(omega_norm[["iter"]])
  norm_names <- unique(betas_norm[["term"]])

  # Prepare all possible combinations
  combos <- expand.grid("Var1" = norm_names, "Var2" = norm_names)
  combos[["combo_raw"]] <- paste(combos[["Var1"]], combos[["Var2"]], sep = "|")
  combos[["combo_unique"]] <- apply(combos, MARGIN = 1L, FUN = function(row) {
    paste(sort(c(row[["Var1"]], row[["Var2"]])), collapse = "|")
  })

  # Match them in the omega df
  omega_norm[["combo_raw"]] <- paste(omega_norm[["Var1"]], omega_norm[["Var2"]], sep = "|")
  omega_norm[["term"]] <- combos[["combo_unique"]][match(
    omega_norm[["combo_raw"]], combos[["combo_raw"]]
  )]
  omega_norm[["combo_raw"]] <- omega_norm[["Var1"]] <- omega_norm[["Var2"]] <- NULL
 # omega_norm <- omega_norm[!duplicated(omega_norm), ]
  #omega_norm <- unique(omega_norm)
  omega_norm <- dplyr::distinct(omega_norm)

  # Return list ready for plotting
  res <- list(
    "betas_subst" = betas_subst,
    "betas_norm" = betas_norm,
    "covar_norm" = omega_norm
  )

  return(res)
}

assess_jomo_chain <- function(MCMCchain_object,
                              formula_subst,
                              dat,
                              intercept = FALSE,
                              ...) {

  param_dfs <- prepare_param_dfs(
    MCMCchain_object = MCMCchain_object,
    formula_subst = formula_subst,
    dat = dat,
    intercept = intercept
  )

  extra_args <- list(...)

  # Now prepare MCMC objects - first subst betas
  long_dats <- lapply(param_dfs, function(df, extra_mcmc_args) {
    # long_df <- reshape(
    #   data = df,
    #   direction = "wide",
    #   idvar = "iter",
    #   timevar = "term",
    #   v.names = "value",
    #   sep = "."
    # )[, -1L] # remove the iter column, not a parameter
    long_df <- tidyr::pivot_wider(df, names_from = "term", values_from = "value")
    long_df$iter <- NULL

    do.call(coda::mcmc, args = c(list("data" = long_df), extra_mcmc_args))
  }, extra_mcmc_args = extra_args)

  return(long_dats)
}

run_example <- FALSE

if (run_example) {

  form <- as.formula(Surv(time, status) ~ measure + sex + I(measure^2))
  dat <- subset(surdata, select = -id)
  imps <- jomo.coxph.MCMCchain(
    formula = form,
    data = dat,
    nburn = 50L,
    #beta.start = matrix(1000, nrow = 1, ncol = 2), # Starting vals don't really work properly...
    #l1cov.start = diag(100, nrow = 2)
  )

  obj <- assess_jomo_chain(imps, form, subset(surdata, select = -id))

  obj$betas_subst |>  head()
  obj$betas_norm |>  head()
  obj$covar_norm |>  head()

  # Traces
  bayesplot::mcmc_trace(x = obj$betas_subst)
  bayesplot::mcmc_trace(x = obj$betas_norm)
  bayesplot::mcmc_trace(x = obj$covar_norm)

  # ACFs
  #lapply(obj, bayesplot::mcmc_acf_bar)
  bayesplot::mcmc_acf_bar(x = obj$betas_subst)
  bayesplot::mcmc_acf_bar(x = obj$betas_norm)
  bayesplot::mcmc_acf_bar(x = obj$covar_norm)

  # Other
  bayesplot::mcmc_dens(x = obj$betas_subst)
  bayesplot::mcmc_intervals(x = obj$betas_subst)
  bayesplot::mcmc_pairs(x = obj$betas_norm)
}
