# Libraries neede (can probably skip naniar)
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

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2))

# Probably only set seed per imputation method - set.seed(...)

# Change later to rds or whatever
dat <- fst::read_fst("data-raw/dat-mds_admin-cens.fst")


# Visualisation -----------------------------------------------------------


#naniar::gg_miss_upset(dat, nsets = 9)
#mice::md.pattern(x = dat, rotate.names = TRUE)
#naniar::vis_miss(x = dat)


# Prepare imputations -----------------------------------------------------


dat$EFS_ind <- as.numeric(dat$ci_s_allo1 > 0)

# Prepare model formula
outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1", "EFS_ind")
predictors <- colnames(dat)[!(colnames(dat) %in% outcomes)]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = predictors)
mod_formula



# CCA ---------------------------------------------------------------------


mod_CCA <- coxph(mod_formula, data = dat)
mod_CCA$coefficients |>  length() # 19 coefficients
mod_CCA$coefficients

# General imputation settings
m <- 3
iters <- 2


future::plan(future::sequential())



# Jomo --------------------------------------------------------------------


# Re-write this part now!!

# Check workflow before continuing https://journal.r-project.org/archive/2019/RJ-2019-034/RJ-2019-034.pdf
# https://www.ucl.ac.uk/population-health-sciences/sites/population_health_sciences/files/quartagno_1.pdf
# https://rmisstastic.netlify.app/tutorials/erler_course_multipleimputation_2018/erler_practical_miadvanced_2018#imputation_using_jomo


# use .MCMCchain to check convergence
# Note jomo used all other vars in df as auxiliary ones, so exclude!!

imps_burn1 <- readRDS("jomo-test-chain.rds")
dim(imps_burn)
jomo.coxph.MCMCchain(beta.start = , l1cov.start =
)
n_burn_iters <- dim(imps_burn$collectbeta)[3]
imps_burn$collectbeta[, , n_burn_iters]
imps_burn$collectomega[, , n_burn_iters]

source("jomo-helpers.R")
imps_burn <- readRDS("jomo-test-chain_start-vals.rds")


obj <- assess_jomo_chain(MCMCchain_object = imps_burn,
                         formula_subst = mod_formula,
                         dat = dat)

bayesplot::mcmc_trace(x = obj$betas_subst)
bayesplot::mcmc_trace(x = obj$betas_norm)
bayesplot::mcmc_pairs(x = obj$betas_norm, regex_pars = "age")
bayesplot::mcmc_trace(
  x = obj$covar_norm,
  regex_pars = "karnof"
)
bayesplot::mcmc_areas_ridges(obj$betas_subst)
bayesplot::mcmc_acf(x = obj$betas_subst, lags = 300)
MCMCvis::MCMCsummary(obj$betas_subst)
test_ggs <- ggmcmc::ggs(obj$betas_subst)
library(ggmcmc)
p <- ggs_running(ggs(obj$covar_norm), family = "cytog")
ggs_traceplot(ggs(obj$covar_norm), family = "cytog")
p + facet_wrap(~ Parameter, ncol = 4) +
  geom_line(size = 1)

ggs_traceplot(ggs(obj$covar_norm), family = "cytog") +
  facet_wrap(~ Parameter, ncol = 4, scales = "free_y")

ggs_autocorrelation(ggs(obj$covar_norm), family = "cytog") +
  facet_wrap(~ Parameter, ncol = 4, scales = "free_y")

ggs_running(ggs(obj$covar_norm), family = "karnof") +
  facet_wrap(~ Parameter, ncol = 4, scales = "free_y") +
  geom_line(size = 1)
# can remove the thinning arg.. only if already thinned..



# JointAI -----------------------------------------------------------------

future::plan(future::sequential)


# This thing is 5 minutes
jointai_imps <- coxph_imp(
  formula = Surv(ci_allo1, EFS_ind) ~ age_allo1_decades,
  n.iter = 5,
  n.adapt = 5,
  data = dat,
  progress.bar = "text"
)

# 10 fookin minutes!!
jointai_imps <- coxph_imp(
  formula = mod_formula,
  n.iter = 100,
  n.adapt = 5,
  data = dat#,
 # progress.bar = "text"
)

jointai_imps$jagsmodel
summary(jointai_imps)
jointai_imps$comp_info$duration
traceplot(jointai_imps, use_ggplot = TRUE)



# MICE --------------------------------------------------------------------


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
mice_imps <- mice(
  m = m, # set globally
  maxit = iters,
  data = dat,
  predictorMatrix = predmat,
  method = meths_mice
)


# SMC-FCS -----------------------------------------------------------------

meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

meths_smcfcs <- smcfcs(
  originaldata = dat,
  smtype = "coxph",
  smformula = deparse1(mod_formula),
  m = 3,
  numit = 3,
  rjlimit = 1000,
  method = meths_smcfcs
)



# Pooling -----------------------------------------------------------------

plot(mice_imps)

EFS_mice_models <- lapply(
  X = mice::complete(mice_imps, action = "all"),
  function(imp) coxph(mod_formula, data = imp)
)

tidy(pool(EFS_mice_models), conf.int = TRUE)
