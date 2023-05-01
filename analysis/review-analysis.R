# Results for systematic review

# Load libraries
library(readxl)
library(tidyverse)
library(hrbrthemes)
library(extrafont)

# Read-in data
extraction_raw <- janitor::clean_names(
  read_xlsx(
    path = "../review-extraction/extraction-sheet.xlsx",
    sheet = 1L,
    skip = 1L
  )
)

# Excluded articles with reasons
excluded_articles <- extraction_raw |>
  filter(str_detect(keep_for_further_extraction, "No*"))

# Reasons for exclusions
table(excluded_articles$if_not_why)

# Subset now just ones used in review
included_articles <- extraction_raw |>
  filter(!str_detect(keep_for_further_extraction, "No*"))

# Make plot of journals
plot_journals <- included_articles |>
  ggplot(aes(reorder(journal, journal, length))) + #fct_infreq also possible
  geom_bar(fill = viridisLite::viridis(1)) +
  stat_count(
    geom = "text",
    colour = "black",
    family = "Roboto Condensed",
    size = 3.5,
    aes(
      label = paste0(
        "n = ", ..count..,
        " (", round(100 * ..count../sum(..count..), 1), "%)"
      ),
      y = ..count.. + 0.5
    ),
    hjust = 0
  ) +
  coord_flip(ylim = c(0, 60)) +
  labs(x = "Journal", y = "Count") +
  theme_ipsum_rc(grid = "X", base_size = 12, axis_title_size = 14)

ggsave(
  plot_journals,
  filename = "figures/journals-overview.svg",
  width = 8,
  height = 6
)

# Exclusions based on missings at population selection
table(included_articles$exclusions_based_on_any_missings)
table(included_articles$exclusions_based_on_any_missings) |>
  prop.table() |>
  round(digits = 2L)

# Types of models - check when they happen at same time
included_articles |>
  separate_rows(multivariable_mv_model_type, sep = "; ") |>
  group_by(multivariable_mv_model_type) |>
  tally() |>
  mutate(prop = round(100 * n / nrow(included_articles), 2))

table(included_articles$multivariable_mv_model_type) |>
  prop.table() |>
  round(digits = 3L)

# Baseline covariates with missings
table(
  included_articles$were_there_baseline_covariates_in_the_mv_model_with_missings,
  useNA = "ifany"
)

table(
  included_articles$were_there_baseline_covariates_in_the_mv_model_with_missings,
  useNA = "ifany"
) |>
  prop.table() |>
  round(digits = 2L)

table(included_articles$were_these_missing_explicity_reported_if_yes_where, useNA = "ifany")
table(included_articles$implicit, included_articles$if_explicit_method_used)

included_articles |>
  separate_rows(were_these_missing_explicity_reported_if_yes_where, sep = "; ") |>
  group_by(were_these_missing_explicity_reported_if_yes_where) |>
  tally()

# Explicit
sum(table(included_articles$if_explicit_method_used))
table(included_articles$if_explicit_method_used, useNA = "ifany")

# Implicitly
table(included_articles$implicit) |> sum()
table(included_articles$implicit, useNA = "ifany")

# Done
included_articles |>
  mutate(software = ifelse(is.na(software), "Unknown", software)) |>
  separate_rows(software, sep = "; ") |>
  group_by(software) |>
  tally()

included_articles |>
  mutate(software = ifelse(is.na(software), "Unknown", software)) |>
  separate_rows(software, sep = "; ") |>
  group_by(software) |>
  tally() |>
  filter(!software %in% c("R", "SAS", "SPSS", "Stata", "Unknown")) |>
  pull(n) |>
  sum()

included_articles |>
  mutate(software = ifelse(is.na(software), "Unknown", software)) |>
  separate_rows(software, sep = "; ") |>
  group_by(title) |>
  tally() |>
  summarise(sum(n >= 2))
