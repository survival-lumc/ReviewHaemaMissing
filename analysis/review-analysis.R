library(readxl)
library(tidyverse)

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

table(excluded_articles$if_not_why)

# Subset now just ones used in review
included_articles <- extraction_raw |>
  filter(!str_detect(keep_for_further_extraction, "No*"))

# Exclusions based on missings at population selection
table(included_articles$exclusions_based_on_any_missings)

# Types of models - check when they happen at same time
included_articles |>
  separate_rows(multivariable_mv_model_type, sep = "; ") |>
  group_by(multivariable_mv_model_type) |>
  tally()


# Baseline covariates with missings
table(
  included_articles$were_there_baseline_covariates_in_the_mv_model_with_missings,
  useNA = "ifany"
)

table(included_articles$were_these_missing_explicity_reported_if_yes_where, useNA = "ifany")

# Explicit
sum(table(included_articles$if_explicit_method_used))
table(included_articles$if_explicit_method_used, useNA = "ifany")

included_articles

# Implicity
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
  group_by(title) |>
  tally() |>
  summarise(sum(n >= 2))

