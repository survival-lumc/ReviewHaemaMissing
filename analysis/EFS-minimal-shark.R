# Libraries
library(jomo)
library(survival)
library(JointAI)
library(mice)
library(smcfcs)
library(future)
library(JointAI)
library(broom)
library(ggplot2)
library(forcats)
library(data.table)
library(tidyverse)
library(patchwork)

#future::availableCores() # with logical = FALSE

cat("Starting script..\n")
old <- Sys.time()
cat(old)

# Set contrasts for ordered factors
options(contrasts = rep("contr.treatment", 2))

# Load-in raw data
dat_raw <- fst::read_fst("data-raw/dat-mds_admin-cens.fst")

# Add EFS indicator and keep only necessary data
dat_raw$EFS_ind <- as.numeric(dat_raw$ci_s_allo1 > 0)
dat <- subset(dat_raw, select = -c(ci_s_allo1, srv_s_allo1, srv_allo1))

# Prepare model formula
predictors <- colnames(dat)[!(colnames(dat) %in% c("ci_allo1", "EFS_ind"))]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = sort(predictors))

dat_missind <- dat
dat_missind[, predictors] <- lapply(dat_missind[, predictors], function(col) {
  if (is.factor(col)) {
    fct_explicit_na(col)
  } else {
    col[is.na(col)] <- 0L
    col
  }
})
dat_missind$`agedonor_allo1_decades(Missing)` <- ifelse(dat_missind$agedonor_allo1_decades == 0, 1, 0)
mod_formula_missind <- update(mod_formula, . ~ . + `agedonor_allo1_decades(Missing)`)
mod_missind <- coxph(mod_formula_missind, data = dat_missind)
summ_missind_raw <- tidy(mod_missind, conf.int = TRUE, exponentiate = TRUE)
summ_missind <- summ_missind_raw[grep(summ_missind_raw$term, pattern = "(Missing)", invert = TRUE), ]

# CCA
#mod_CCA <- coxph(mod_formula, data = dat)


# FCS approaches ----------------------------------------------------------

cat("\nStarting mice..\n")
print(Sys.time() - old)

# General settings
m <- 20L
cycles <- 15L

# Prepare imputation models
dat$cumhaz <- nelsonaalen(data = dat, timevar = ci_allo1, statusvar = EFS_ind)
predmat <- make.predictorMatrix(dat)
predmat[] <- 0L
predmat[, c(predictors, "EFS_ind", "cumhaz")] <- 1L
diag(predmat) <- 0L

# Prepare methods - replace PMM?
meths_mice <- make.method(dat)

imps_mice <- parlmice(
  #m = m, # set globally
  n.imp.core = 5L,
  maxit = cycles,
  n.core = 4L, # adjust this, number of imps is n.core * n.imp.core
  data = dat,
  predictorMatrix = predmat,
  method = meths_mice,
  seed = 97868
)

cat("\nStarting smcfcs..\n")
print(Sys.time() - old)

meths_smcfcs <- make.method(dat, defaultMethod = c("norm", "logreg", "mlogit", "podds"))

imps_smcfcs <- smcfcs.parallel(
  originaldata = dat,
  smtype = "coxph",
  smformula = deparse1(mod_formula),
  m = m,
  numit = cycles,
  n_cores = 4L,
  rjlimit = 5000,
  method = meths_smcfcs,
  seed = 45684
)


# JM approaches -----------------------------------------------------------


# Set seed now for remaining methods
set.seed(3002764)


cat("\nStarting jomo..\n")
print(Sys.time() - old)

# Keep only variables in substantive model
dat_jomo <- subset(dat, select = -cumhaz)

# Also crashes when out.iter = nbetween?
imps_jomo <- jomo.coxph(
  formula = mod_formula,
  data = dat_jomo,
  nburn = 5000L,
  nbetween = 250L,
  out.iter = 100L,
  nimp = m
)

cat("\nStarting jointai..\n")
print(Sys.time() - old)

plan(multicore, workers = 4)

mod <- coxph_imp(
  mod_formula,
  n.chains = 4,
  data = dat,
  n.iter = 2000,
  n.adapt = 200,
  monitor_params = list(analysis_main = TRUE, alphas = TRUE)
)

mod$comp_info # time
print(Sys.time() - old)

imps_all <- list(
  "mice_obj" = imps_mice,
  "smcfcs_obj" = imps_smcfcs,
  "jomo_obj" = imps_jomo,
  "jointai_obj" = mod
)

saveRDS(imps_all, file = "imps_all.rds")

cat("\nEnd!")

# Can lower df_basehaz later, too flexible atm



# Analysis ----------------------------------------------------------------


mod_CCA <- coxph(mod_formula, data = dat)
imps_all <- readRDS("data/imps_all_2.rds")

EFS_mice_models <- lapply(
  X = mice::complete(imps_all$mice_obj, action = "all"),
  function(imp) coxph(mod_formula, data = imp)
)

EFS_smcfcs_models <- lapply(
  X = imps_all$smcfcs_obj$impDatasets,
  function(imp) coxph(mod_formula, data = imp)
)

jomo_imps <- subset(imps_all$jomo_obj, subset = Imputation != 0)

EFS_jomo_models <- lapply(
  X = split(jomo_imps, f = jomo_imps$Imputation),
  function(imp) coxph(mod_formula, data = imp)
)

#EFS_CCA <- tidy(mod_CCA)
traceplot(imps_all$jointai_obj,
          subset = c(analysis_main = FALSE, other_models = TRUE))
          #subset = c(analysis_main = FALSE, other_models = TRUE))

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

summ_all <- dplyr::bind_rows(
  tidy_jointai,
  cbind(tidy(pool(EFS_mice_models), conf.int = TRUE), "method" = "MICE"),
  cbind(tidy(pool(EFS_smcfcs_models), conf.int = TRUE), "method" = "SMC-FCS"),
  #cbind(tidy(pool(EFS_jomo_models), conf.int = TRUE), "method" = "jomo"),
  cbind(summ_missind, "method" = "Missing\nindicator"),
  cbind(tidy(mod_CCA, conf.int = TRUE), "method" = "CCA")
)

library(dplyr)
summ_all |>
  select(term, method, estimate, conf.low, conf.high) |>
  ggplot(aes(term, estimate, group = method)) +
  geom_linerange(
    aes(
      ymin = conf.low,
      ymax = conf.high,
      xmin = term,
      xmax = term,
      col = method
    ),
    position = position_dodge(width = 0.75),
    size = 0.5,
    na.rm = TRUE
  ) +
  geom_point(
    aes(col = method, shape = method),
    position = position_dodge(width = 0.75),
    size = 1.75, #1.25,
    na.rm = TRUE
  ) +
  coord_flip() +
  theme_minimal()


library(ggforestplot)
summ_all$method

p <- summ_all |>
  mutate(
    method = factor(
      method,
      levels = rev(c("CCA", "Missing\nindicator", "MICE", "SMC-FCS", "JointAI"))
    )
  ) |>
  forestplot(
    name = term,
    logodds = TRUE,
    estimate = estimate,
    se = std.error, # not true for bayesian!
    colour = method,
    shape = method,
    xlab = "Hazard ratio (95% CI)",
    xtickbreaks = c(0.75, 1, 1.5, 2)
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(col = "Method", shape = "Method") +
  theme_forest(base_size = 20)

p

png(p, filename = "forest-raw.png", width = 12, height = 12, res = 300, units = "in")
p
dev.off()


library("forestploter")
dt <- read.csv(system.file("extdata", "example_data.csv", package = "forestploter"))


merge(summ_all, dictionary_df) |> View()





# Test with old function --------------------------------------------------

#remotes::install_github("survival-lumc/CauseSpecCovarMI", build_vignettes = FALSE)

library(CauseSpecCovarMI)
#load(system.file("R", "data_dictionary.rda", package = "CauseSpecCovarMI"))
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
  summ$`97.5 %` <- summ$conf.high
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
  filename = "./forest_EFS.pdf",
  plot = p,
  width = 10,
  height = 13, dpi = 300
)

ggplot2::ggsave(
  filename = "./forest_EFS.png",
  plot = p,
  width = 10,
  height = 13, dpi = 300
)

