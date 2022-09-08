# Seed
set.seed(7367)
source("jomo-helpers.R")

# Read-in
dat <- data.frame(fst::read_fst("data-raw/dat-mds_admin-cens.fst"))

# Add indicator
dat$EFS_ind <- as.numeric(dat$ci_s_allo1 > 0)
dat$karnofsk_allo1 <- factor(dat$karnofsk_allo1, ordered = FALSE)
dat$cytog_threecat <- factor(dat$cytog_threecat, ordered = FALSE)
dat$hctci_risk <- factor(dat$hctci_risk, ordered = FALSE)

# Prepare formula
outcomes <- c("ci_s_allo1", "ci_allo1", "srv_s_allo1", "srv_allo1", "EFS_ind")
predictors <- colnames(dat)[!(colnames(dat) %in% outcomes)]
mod_formula <- reformulate(response = "Surv(ci_allo1, EFS_ind)", termlabels = predictors)
mod_formula

# Subset
dat_sub <- subset(dat, select = -c(ci_s_allo1, srv_allo1, srv_s_allo1))

jomo_imps_chain <- jomo.coxph.MCMCchain(
  formula = mod_formula,
  data = dat_sub,
  nburn = 100,
  out.iter = 10
)


# Specific priors ---------------------------------------------------------

#
#initial_chains <- readRDS("jomo-test-chain.rds")


# Attempt for specific priors/starting values.. :
# n_iter <- dim(initial_chains$collectbeta)[3]
# df <- dim(initial_chains$collectbeta)[2]
# beta_init <- initial_chains$collectbeta[, , n_iter]
# covmat_init <- initial_chains$collectomega[, , n_iter]
# covmat_scaled<- covmat_init * df

