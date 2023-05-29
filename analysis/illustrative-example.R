# Handling missing covariate data in clinical studies in haematology (2023)
# E. F. Bonneville et al.

# Code for illustrative example (section 5)


# Packages and general settings -------------------------------------------


# Packages needed:
library(survival) # For fitting Cox models in each imputed dataset
library(mice) # MICE approach
library(smcfcs) # SMC-FCS approach
library(JointAI) # Fully Bayesian approach
library(future) # In order to use parallel computation in the Bayesian approach
library(ggplot2) # Plotting
library(data.table) # Data manipulation
library(forcats) # Handling factors for the missing indicator method
library(patchwork) # Extra plotting
library(broom) # Extracting model summaries

## General settings:

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2))

# Detect number of available cores for parallel computation
# (analysis reproducible with num_cores = 4)
num_cores <- parallelly::availableCores(logical = FALSE) - 1L


# Read-in data ------------------------------------------------------------


# Load-in raw data
dat_raw <- fst::read_fst("data-raw/dat-mds_admin-cens.fst")

# Add EFS indicator and keep only necessary data
dat_raw$EFS_ind <- as.numeric(dat_raw$ci_s_allo1 > 0)
dat <- subset(dat_raw, select = -c(ci_s_allo1, srv_s_allo1, srv_allo1))

# Check remaining columns names
colnames(dat)
# [1] "ci_allo1"               "match_allo1_1"          "mdsclass"               "donorrel"
# [5] "karnofsk_allo1"         "crnocr"                 "cmv_combi_allo1_1"      "cytog_threecat"
# [9] "hctci_risk"             "agedonor_allo1_decades" "age_allo1_decades"      "EFS_ind"

# Prepare analysis model formula
predictors <- colnames(dat)[!(colnames(dat) %in% c("ci_allo1", "EFS_ind"))]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = sort(predictors))
mod_formula
# Surv(ci_allo1, EFS_ind) ~ age_allo1_decades + agedonor_allo1_decades +
#   cmv_combi_allo1_1 + crnocr + cytog_threecat + donorrel +
#   hctci_risk + karnofsk_allo1 + match_allo1_1 + mdsclass


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
print(summ_missind_raw)
# A tibble: 26 × 7
# term                                   estimate std.error statistic  p.value conf.low conf.high
# <chr>                                     <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#   1 age_allo1_decades                         1.11     0.0156     6.55  5.62e-11    1.07       1.14
# 2 agedonor_allo1_decades                    1.07     0.0188     3.61  3.07e- 4    1.03       1.11
# 3 cmv_combi_allo1_1-/+                      0.978    0.0671    -0.332 7.40e- 1    0.857      1.12
# 4 cmv_combi_allo1_1+/-                      1.12     0.0505     2.28  2.25e- 2    1.02       1.24
# 5 cmv_combi_allo1_1+/+                      1.06     0.0464     1.28  2.00e- 1    0.969      1.16
# 6 cmv_combi_allo1_1(Missing)                1.09     0.0636     1.38  1.66e- 1    0.964      1.24
# 7 crnocrnoCR                                1.33     0.0404     6.96  3.29e-12    1.22       1.43
# 8 crnocrUntreated/not aimed at remission    1.04     0.0469     0.875 3.81e- 1    0.950      1.14
# 9 crnocr(Missing)                           1.30     0.0942     2.79  5.34e- 3    1.08       1.56
# 10 cytog_threecatpoor(4)                     1.33     0.0793     3.63  2.82e- 4    1.14       1.56
# # … with 16 more rows

# Omit the missing indicator coefficients from output
summ_missind <- summ_missind_raw[grep(summ_missind_raw$term, pattern = "(Missing)", invert = TRUE), ]


# MI approaches -----------------------------------------------------------


# General settings
m <- 100L # Number of imputed datasets
cycles <- 15L # Number of iterations

# -- MICE

## Prepare imputation models

# Add cumulative hazard to dataset
dat$cumhaz <- nelsonaalen(data = dat, timevar = ci_allo1, statusvar = EFS_ind)

# Specify what variables to be included inside each imputation model
predmat <- make.predictorMatrix(dat)
predmat[] <- 0L
predmat[, c(
  predictors, # Predictors from the analysis model
  "EFS_ind", # Event indicator
  "cumhaz" # Cumulative hazard
)] <- 1L
diag(predmat) <- 0L # Variable is not included in its own imputation model

# Prepare imputation methods
# (Note predictive mean matching, PMM, is used for continuous variables)
# Check ?mice for available imputation methods..
meths_mice <- make.method(dat)
meths_mice
# ci_allo1          match_allo1_1               mdsclass               donorrel
# ""              "polyreg"                     ""                     ""
# karnofsk_allo1                 crnocr      cmv_combi_allo1_1         cytog_threecat
# "polr"              "polyreg"              "polyreg"                 "polr"
# hctci_risk agedonor_allo1_decades      age_allo1_decades                EFS_ind
# "polr"                  "pmm"                     ""                     ""
# cumhaz
# ""

# Run MICE using parallel version of mice() function
imps_mice <- parlmice(
  n.imp.core = ceiling(m / num_cores), # Number of imputations per core
  maxit = cycles,
  n.core = num_cores,
  data = dat,
  predictorMatrix = predmat,
  method = meths_mice,
  cluster.seed = 97868 # For reproducibility
)

# -- SMC-FCS


# SMC-FCS uses slightly difference codes for imputation methods:
meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

# Run SMC-FCS imputations (no need for predictorMatrix)
imps_smcfcs <- smcfcs.parallel(
  originaldata = dat,
  smtype = "coxph",
  smformula = deparse1(mod_formula),
  m = m,
  numit = cycles,
  n_cores = num_cores,
  rjlimit = 10000L, # Set higher if warnings occur
  method = meths_smcfcs,
  seed = 45684
)


# JM approaches -----------------------------------------------------------


# Set seed now for remaining method
set.seed(3002764)

# Set plan to multisession if on windows
plan(multicore, workers = num_cores)

# Run Fully Bayesian imputations
imps_jointai <- coxph_imp(
  mod_formula,
  n.chains = 4L,
  data = dat,
  n.iter = 2000L,
  n.adapt = 200L,
  monitor_params = list(analysis_main = TRUE, alphas = TRUE)
)


# Put all imputation objects together in a list
imps_all <- list(
  "mice_obj" = imps_mice,
  "smcfcs_obj" = imps_smcfcs,
  "jointai_obj" = imps_jointai
)

# Save in order to avoid re-running imputations
saveRDS(imps_all, file = "./data/imps_all.rds")


# Complete-case analysis --------------------------------------------------


mod_CCA <- coxph(mod_formula, data = dat)


# Pooling models ----------------------------------------------------------


# Load saved objects
#imps_all <- readRDS("data/imps_all.rds")

# Fit models in imputed datasets for both MICE and SMC-FCS, using
# previously created substantive model formula
EFS_mice_models <- lapply(
  X = mice::complete(imps_all$mice_obj, action = "all"),
  function(imp) coxph(mod_formula, data = imp)
)

EFS_smcfcs_models <- lapply(
  X = imps_all$smcfcs_obj$impDatasets,
  function(imp) coxph(mod_formula, data = imp)
)

# Check convergence of imputations too:
#traceplot(
#  imps_all$jointai_obj,
#  subset = c(analysis_main = FALSE, other_models = TRUE)
#)
#plot(imps_all$mice_obj)
#plot(imps_all$smcfcs_obj)


# Summarise results -------------------------------------------------------


# (This part reproduces Figure in manuscript)

# Summarise JointAI object
EFS_jointai <- summary(imps_all$jointai_obj)$res$`Surv(ci_allo1, EFS_ind)`$regcoef
tidy_jointai <- cbind.data.frame(
  "term" = rownames(EFS_jointai),
  "estimate" = exp(EFS_jointai[, "Mean"]),
  "std.error" = EFS_jointai[, "SD"],
  "p.value" = EFS_jointai[, "tail-prob."],
  "conf.low" = exp(EFS_jointai[, "2.5%"]),
  "conf.high" = exp(EFS_jointai[, "97.5%"])
)

# Helper functions/objects for the plotting
load("./data-raw/data_dictionary.rda")
source("./R/forest-helper.R")

# Combine results
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

# Make forest plot
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

# Save
ggplot2::ggsave(
  filename = "./figures/forest_EFS.pdf",
  plot = p,
  width = 10,
  height = 13, dpi = 300
)
