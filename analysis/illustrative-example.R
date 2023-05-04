# Libraries
library(survival)
library(JointAI)
library(mice)
library(smcfcs)
library(future)
library(ggplot2)
library(forcats)
library(data.table)
library(patchwork)
library(broom)
#library(tidyverse)

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2))

# Detect number of available cores for the analysis
num_cores <- parallelly::availableCores(logical = FALSE)


# Read-in data ------------------------------------------------------------


# Load-in raw data
dat_raw <- fst::read_fst("data-raw/dat-mds_admin-cens.fst")

# Add EFS indicator and keep only necessary data
dat_raw$EFS_ind <- as.numeric(dat_raw$ci_s_allo1 > 0)
dat <- subset(dat_raw, select = -c(ci_s_allo1, srv_s_allo1, srv_allo1))

# Prepare model formula
predictors <- colnames(dat)[!(colnames(dat) %in% c("ci_allo1", "EFS_ind"))]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = sort(predictors))


# Missing indicator method (MIM) ------------------------------------------


# Make a dataset specifically for missing indicator method (MIM)
dat_missind <- dat

# Create the missing indicators (note all character columns should already
# be set to factors)
dat_missind[, predictors] <- lapply(dat_missind[, predictors], function(col) {
  if (!is.factor(col)) {
    col[is.na(col)] <- 0L; col # Set NAs to zero for continuous variables
  } else fct_explicit_na(col) # Add extra level for categorical variables
})

# Make indicator for donor age (only continuous variable with missings)
dat_missind$`agedonor_allo1_decades(Missing)` <- ifelse(dat_missind$agedonor_allo1_decades == 0, 1, 0)

# Update formula to include indicators
mod_formula_missind <- update(mod_formula, . ~ . + `agedonor_allo1_decades(Missing)`)

# Run MIM, and summarise
mod_missind <- coxph(mod_formula_missind, data = dat_missind)
summ_missind_raw <- tidy(mod_missind, conf.int = TRUE, exponentiate = TRUE)

# Omit the missing indicator coefficients from output
summ_missind <- summ_missind_raw[grep(summ_missind_raw$term, pattern = "(Missing)", invert = TRUE), ]


# MI approaches -----------------------------------------------------------


# General settings
m <- 10 #100L
cycles <- 50L

# -- MICE

# Prepare imputation models (add cumulative hazard to dataset)
dat$cumhaz <- nelsonaalen(data = dat, timevar = ci_allo1, statusvar = EFS_ind)
predmat <- make.predictorMatrix(dat)
predmat[] <- 0L
predmat[, c(predictors, "EFS_ind", "cumhaz")] <- 1L
diag(predmat) <- 0L
predmat

# Prepare methods
meths_mice <- make.method(dat)
meths_mice

imps_mice <- parlmice(
  n.imp.core = ceiling(m / num_cores),
  maxit = cycles,
  n.core = num_cores,
  data = dat,
  predictorMatrix = predmat,
  method = meths_mice,
  seed = 97868
)

# -- SMC-FCS

meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

imps_smcfcs <- smcfcs.parallel(
  originaldata = dat,
  smtype = "coxph",
  smformula = deparse1(mod_formula),
  m = m,
  numit = cycles,
  n_cores = num_cores,
  rjlimit = 10000L,
  method = meths_smcfcs,
  seed = 45684
)


# JM approaches -----------------------------------------------------------


# Set seed now for remaining methods
set.seed(3002764)

# Set plan to multisession if on windows..
plan(multicore, workers = num_cores)

mod <- coxph_imp(
  mod_formula,
  n.chains = 4L,
  data = dat,
  n.iter = 200L, #2000L,
  n.adapt = 200L,
  monitor_params = list(analysis_main = TRUE, alphas = TRUE)
)


# Put all imputation objects together
imps_all <- list(
  "mice_obj" = imps_mice,
  "smcfcs_obj" = imps_smcfcs,
  "jointai_obj" = mod
)

saveRDS(imps_all, file = "data/imps_all.rds")

stop()


# Complete-case analysis --------------------------------------------------


mod_CCA <- coxph(mod_formula, data = dat)


# Pooling models ----------------------------------------------------------


# Load saved objects
#imps_all <- readRDS("data/imps_all_2.rds")

# Fit models in imputed datasets (FIX: mice only 20 impdats!!)
EFS_mice_models <- lapply(
  X = mice::complete(imps_all$mice_obj, action = "all"),
  function(imp) coxph(mod_formula, data = imp)
)

EFS_smcfcs_models <- lapply(
  X = imps_all$smcfcs_obj$impDatasets,
  function(imp) coxph(mod_formula, data = imp)
)

# JointAI convergence good!
traceplot(
  imps_all$jointai_obj,
  subset = c(analysis_main = FALSE, other_models = TRUE)
)


# Summarise results -------------------------------------------------------


EFS_jointai <- summary(imps_all$jointai_obj)$res$`Surv(ci_allo1, EFS_ind)`$regcoef
tidy_jointai <- cbind.data.frame(
  "term" = rownames(EFS_jointai),
  "estimate" = exp(EFS_jointai[, "Mean"]),
  "std.error" = EFS_jointai[, "SD"],
  "p.value" = EFS_jointai[, "tail-prob."],
  "conf.low" = exp(EFS_jointai[, "2.5%"]),
  "conf.high" = exp(EFS_jointai[, "97.5%"])#,
  #"method" = "JointAI"
)

load("data-raw/data_dictionary.rda")
source("R/forest-helper.R")

results_raw <- list(
  "JointAI" = tidy_jointai,
  "MICE" = tidy(pool(EFS_mice_models), conf.int = TRUE, exponentiate = TRUE),
  "SMC-FCS" = tidy(pool(EFS_smcfcs_models), conf.int = TRUE, exponentiate = TRUE),
  "Missing\nindicator" = summ_missind,
  "CCA" = tidy(mod_CCA, conf.int = TRUE, exponentiate = TRUE)
)

results <- lapply(results_raw, function(summ) {
  summ$`2.5 %` <- summ$conf.low;
  summ$`97.5 %` <- summ$conf.high;
  summ
})

p <- ggplot_grouped_forest(
  dat = data.table(dat),
  dictionary = dictionary_df,
  results = results,
  event = "Event-free survival",
  form = mod_formula,
  breaks_x = c(0.5, 1, 1.5, 2),
  lims_x = c(0.65, 2.5)
) +
  ggplot2::theme(
    legend.margin = ggplot2::margin(t = 0, r = 0, b = -1.5, l = 0, unit = "cm")
  )

ggplot2::ggsave(
  filename = "./figures/forest_EFS.pdf",
  plot = p,
  width = 10,
  height = 13, dpi = 300
)

ggplot2::ggsave(
  filename = "./figures/forest_EFS.png",
  plot = p,
  width = 10,
  height = 13, dpi = 300
)
