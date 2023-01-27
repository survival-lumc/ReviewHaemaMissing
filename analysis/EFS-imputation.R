# Libraries needed (can probably skip naniar)
library(naniar)
library(broom)
library(survival)
library(mice)
library(smcfcs)
library(jomo)
library(JointAI)
library(coda)
library(bayesplot)
library(ggplot2)
library(future)

# Source helpers
source("R/jomo-helpers.R")

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2))

# Probably only set seed per imputation method - set.seed(...)

# Change later to rds or whatever - check also whether we want admin cens one
dat_raw <- fst::read_fst("data-raw/dat-mds_admin-cens.fst")


# Prepare data ------------------------------------------------------------


# Note all ordered/unordered already set at this point..

# Make EFS indicator, and subset to keep only relevant variables
# .. potentially bring back later if cuminc plots needed
dat_raw$EFS_ind <- as.numeric(dat_raw$ci_s_allo1 > 0)
dat <- subset(dat_raw, select = -c(ci_s_allo1, srv_s_allo1, srv_allo1))

# Prepare model formula
predictors <- colnames(dat)[!(colnames(dat) %in% c("ci_allo1", "EFS_ind"))]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = sort(predictors)) # maybs unsort later
mod_formula


# Visualisation -----------------------------------------------------------


#naniar::gg_miss_upset(dat, nsets = 9)
#mice::md.pattern(x = dat, rotate.names = TRUE)
#naniar::vis_miss(x = dat)


# CCA ---------------------------------------------------------------------


mod_CCA <- coxph(mod_formula, data = dat)
tidy(mod_CCA)


# FCS approach 1: MICE ----------------------------------------------------


# General settings
m <- 2L
cycles <- 2L

# Prepare imputation models
dat$cumhaz <- nelsonaalen(data = dat, timevar = ci_allo1, statusvar = EFS_ind)
predmat <- make.predictorMatrix(dat)
predmat[] <- 0L
predmat[, c(predictors, "EFS_ind", "cumhaz")] <- 1L
diag(predmat) <- 0L
predmat

# Prepare methods
meths_mice <- make.method(dat)
meths_mice # replace the pmm?

# use parlmice/try future version?
imps_mice <- mice(
  m = m, # set globally
  maxit = cycles,
  data = dat,
  predictorMatrix = predmat,
  method = meths_mice
)

plot(imps_mice)


# FCS approach 2: SMC-FCS -------------------------------------------------


meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

imps_smcfcs <- smcfcs(
  originaldata = dat,
  smtype = "coxph",
  smformula = deparse1(mod_formula),
  m = m,
  numit = cycles,
  rjlimit = 5000, # Edit probably
  method = meths_smcfcs
)

plot(imps_smcfcs)


# Joint model approaches: jomo --------------------------------------------


# Keep only variables in substantive model
dat_jomo <- subset(dat, select = -cumhaz)

# First dry run to check for convergence
inits_jomo <- jomo.coxph.MCMCchain(
  formula = mod_formula,
  data = dat_jomo,
  nburn = 10L,
  out.iter = 10L # adjust these
)

# Convergence code here..
nburn_samples <- assess_jomo_chain(
  MCMCchain_object = inits_jomo,
  formula_subst = mod_formula,
  dat = dat_jomo
)

bayesplot::mcmc_trace(x = obj$betas_subst)
bayesplot::mcmc_trace(x = obj$betas_norm)
bayesplot::mcmc_pairs(x = obj$betas_norm, regex_pars = "age")
MCMCvis::MCMCsummary(obj$betas_subst)
ggmcmc::ggs_running(ggmcmc::ggs(obj$covar_norm), family = "cytog")

# Actual imputations
nburn_init <- dim(inits_jomo$collectbeta)[3]
imps_jomo <- jomo.coxph(
  formula = mod_formula,
  data = dat_jomo,
  beta.start = matrix(inits_jomo$collectbeta[, , nburn_init], nrow = 1L),
  l1cov.start = inits_jomo$collectomega[, , nburn_init],
  nburn = 5000L, # something smaller than the above
  nbetween = 500L,
  out.iter = 100L,
  nimp = m
)


# Joint model approaches: JointAI -----------------------------------------


future::plan(sequential) # change to multicore

# 10 fookin minutes!!
jointai_imps <- coxph_imp(
  formula = mod_formula,
  n.iter = 100,
  n.adapt = 5,
  data = dat
)

#jointai_imps$jagsmodel
summary(jointai_imps)
jointai_imps$comp_info$duration
traceplot(jointai_imps, use_ggplot = TRUE)


# Pooling -----------------------------------------------------------------


EFS_mice_models <- lapply(
  X = mice::complete(mice_imps, action = "all"),
  function(imp) coxph(mod_formula, data = imp)
)

tidy(pool(EFS_mice_models), conf.int = TRUE)
