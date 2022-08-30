library(readxl)
library(tidyverse)

# Read-in sheet
dat_journos <- readxl::read_xlsx(
  path = "../EdouardBonneville_JCR_JournalResults_05_2022.xlsx"
)

# Keep only relevant journals
chosen_journos <- dat_journos |>
  filter(Include == "Yes" & !is.na(Include))

name_query <- paste(paste0("'", chosen_journos$`Journal name`, "'"), collapse = " or ")
test <- gsub(name_query, pattern = "'", replacement = '"')
print(test, quote = FALSE)

# Try also issns
issn_query <- paste(paste0("'", chosen_journos$ISSN, "'"), collapse = " OR ")
issn_test <- gsub(issn_query, pattern = "'", replacement = '"')
print(issn_test, quote = FALSE)
cat(paste(chosen_journos$ISSN, collapse = " OR "))

# Abbrevations are just clarivate abbrevs! use journal names!!
JIF_geq_two <- paste(paste0("'", chosen_journos$`JCR Abbreviation`, "'"), collapse = " or ")
test <- gsub(JIF_geq_two, pattern = "'", replacement = '"')

# Web of science format
paste(paste0("'", chosen_journos$`Journal name`, "'"), collapse = " OR ")

# Try for scholar
paste(paste0("source:'", chosen_journos$`JCR Abbreviation`, "'"), collapse = " OR ")

# Add cox (imputation|missing) just prior
clipr::write_clip(paste(paste0("source:'", chosen_journos$`JCR Abbreviation`, "'"), collapse = " OR "))



# Try importing titles ----------------------------------------------------

test_corpus <- readxl::read_xls("../test-ovid-export.xls", skip = 1)
print(
  paste(paste0('"', test_corpus$TI[1:2], '"'), collapse = " OR "),
  quote = FALSE
)

clipr::write_clip(
  paste(
    paste0("intitle:'", gsub(test_corpus$TI[1:3], pattern = "\\.", replacement = ""), "'"),
    collapse = " OR "
  )
)


# Try other things
paste(
  paste0("intitle:'", gsub(test_corpus$TI[1:2], pattern = "\\.", replacement = ""), "'"),
    collapse = " OR "
)

paste(
  paste0("author:", stringr::word(test_corpus$AU, start = 1)),
  collapse = " OR "
)

