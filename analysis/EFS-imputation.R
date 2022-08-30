# This will be only shared file


library(naniar)
library(broom)
library(survival)
library(mice)
library(smcfcs)
library(jomo)
library(JointAI)
#library(Hmisc) aregImpute
#library(Amelia)
# Latest bartlett reference based MI?

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2))

set.seed(79879)

# https://academic.oup.com/aje/article/179/6/764/107562
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-10-112

# Change later to rds or whatever
dat <- fst::read_fst("data-raw/dat-mds_admin-cens.fst")


# Visualisation -----------------------------------------------------------


gg_miss_upset(dat)
mice::md.pattern(x = dat, rotate.names = TRUE)
naniar::vis_miss(x = dat)


# Prepare imputations -----------------------------------------------------


dat$EFS_ind <- as.numeric(dat$ci_s_allo1 > 0)

# Prepare model formula
outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1", "EFS_ind")
predictors <- colnames(dat)[!(colnames(dat) %in% outcomes)]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = predictors)
mod_formula

# General imputation settings
m <- 3
iters <- 2


future::plan(future::sequential())

# Check workflow before continuing https://journal.r-project.org/archive/2019/RJ-2019-034/RJ-2019-034.pdf
# https://www.ucl.ac.uk/population-health-sciences/sites/population_health_sciences/files/quartagno_1.pdf

# use .MCMCchain to check convergence

jomo_imps <- jomo.coxph(
  formula = mod_formula,
  data = dat,
  nimp = 5,
  nburn = 100,
  nbetween = 100
)

library(doFuture)
doFuture::registerDoFuture()
future::plan(future::sequential)

# This thing is 5 minutes
jointai_imps <- coxph_imp(
  formula = Surv(ci_allo1, EFS_ind) ~ age_allo1_decades,
  n.iter = 5,
  n.adapt = 5,
  data = dat,
  progress.bar = "gui"
)

# 10 fookin minutes!!
jointai_imps <- coxph_imp(
  formula = mod_formula,
  n.iter = 1,
  n.adapt = 0,
  data = dat,
  progress.bar = "gui"
)

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
