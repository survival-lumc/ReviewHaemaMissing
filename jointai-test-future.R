library(JointAI)
library(future)


mod1 <- lm_imp(y ~ C1 + C2 + M1 + B1, data = wideDF, n.iter = 10000)

plan(multisession, workers = 3)


system.time({
  plan(sequential)
  mod1 <- lm_imp(
    y ~ C1 + C2 + M1 + B1, data = wideDF,
    n.iter = 100000,
    n.chains = 3
  )
})

system.time({
  plan(multisession, workers = 3)
  mod1 <- lm_imp(
    y ~ C1 + C2 + M1 + B1, data = wideDF,
    n.iter = 100000,
    n.chains = 3
  )
})

microbenchmark::microbenchmark(
  seq = {
    plan(sequential)
    mod1 <- lm_imp(
      y ~ C1 + C2 + M1 + B1, data = wideDF,
      n.iter = 100000,
      n.chains = 3
    )
  },
  threecore = {
    plan(multisession, workers = 3)
    mod1 <- lm_imp(
      y ~ C1 + C2 + M1 + B1, data = wideDF,
      n.iter = 100000,
      n.chains = 3
    )
  },
  times = 5L
)


summary(jointai_obj)
library(JointAI)
library(ggplot2)
jointai_obj
traceplot(
  jointai_obj, ncol = 4, use_ggplot = TRUE,
  subset = list(alphas = FALSE)
) +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_brewer(palette = 'Dark2')

traceplot(
  jointai_obj, ncol = 4, use_ggplot = TRUE,
  subset = list(analysis_main = FALSE, alphas = TRUE)
) +
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_brewer(palette = 'Dark2')

GR_crit(jointai_obj, subset = list(analysis_main = FALSE, alphas = TRUE))

densplot(
  jointai_obj, ncol = 4, use_ggplot = TRUE,
  subset = list(alphas = FALSE)
)

coef(jointai_obj)
