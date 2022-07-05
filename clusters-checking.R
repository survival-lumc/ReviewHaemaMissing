library(readxl)
library(tidyverse)

# Read-in sheet
corpus_raw <- readxl::read_xls(path = "../2022-07-04_export.xls", skip = 1L)

#
corpus <- corpus_raw |>
  mutate(
    AU = str_replace_all(AU, pattern = "(\n\n)", replacement = ","),
    KW = str_replace_all(KW, pattern = "(\n\n)", replacement = ",")
  ) |>
  select(-all_of(c("VN", "FTURL", "XL")))
 # View()

# Basic descriptives
table(corpus$JA)

# Liesbeth
corpus |> filter(str_detect(AU, "Wreede")) # remove
corpus |> filter(str_detect(AU, "Wreede")) |> View()
corpus |> filter(str_detect(AU, "Schetelig"))
corpus |> filter(str_detect(AU, "Polverelli"))


# Checking author clustering
corpus |>
  mutate(n_total_articles = n()) |>
  separate_rows(AU, sep = ",") |>
  group_by(AU, n_total_articles) |>
  tally() |>
  ungroup() |>
  mutate(prop = paste0(round(100 * n / n_total_articles, 2), "%")) |>
  filter(n >= 5) |>
  arrange(desc(n)) |>
  View()

# Check clustering of CIBMTR/EBMT -> check in abstract too
corpus |> filter(str_detect(TI, "EBMT") | str_detect(AB, "EBMT"))
corpus |> filter(str_detect(TI, "CIBMTR") | str_detect(AB, "EBMT"))


# Try capturing malignancy from title -------------------------------------

corpus |> filter(str_detect(TI, "CLL"))
corpus |> filter(str_detect(TI, "AML"))
corpus |> filter(str_detect(TI, "MDS"))
corpus |> filter(str_detect(TI, "HL"))
corpus |> filter(str_detect(TI, "NHL"))
corpus |> filter(str_detect(TI, "MPN"))
corpus |> filter(str_detect(TI, "MM"))
corpus |> filter(str_detect(TI, "ALL"))

# ss ----------------------------------------------------------------------


# Abstract works better than title
corpus_malign <- corpus |>
  mutate(
    malignancy = case_when(
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "CLL|([Cc]hronic [Ll]ymphocytic [Ll]eukemia)")) ~ "CLL",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "CML|([Cc]hronic [My]yeloid [Ll]eukemia)")) ~ "CML",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "AML|([Aa]cute [Mm]yeloid [Ll]eukemia)|([Aa]cute [Mm]yelogenous [Ll]eukemia)")) ~ "AML",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "ALL|([Aa]cute [Ll]ymphocytic [Ll]eukemia)|([Aa]cute [Ll]ymphoblastic [Ll]eukemia)")) ~ "ALL",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "MPN|([Mm]yeloproliferative [Nn]eoplasm*)")) ~ "MPN",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "MM|([Mm]ultiple [Mm]yeloma*)|[Mm]yeloma*")) ~ "Myelomas (incl MM)",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "HL|Hodgkin*")) ~ "HL",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "NHL|non-Hodgkin*")) ~ "NHL",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "[Ll]ymphoma*")) ~ "Lymphomas",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "MDS|([Mm]yelodysplastic [Ss]yndrome*)")) ~ "MDS",
      if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "CLL|([Cc]hronic [Ll]ymphocytic [Ll]eukemia)")) ~ "CLL",
      TRUE ~ "Unclear"
    )
  ) #|>
 # View()
  #group_by(malignancy) |>
  #tally()
  #filter(malignancy == "Unclear") |>
  #View()
str_extract_all(
  string = corpus$TI,
  pattern = "CLL|([Cc]hronic [Ll]ymphocytic [Ll]eukemia)|[Ll]ymphoma*|CLL|([Cc]hronic [Ll]ymphocytic [Ll]eukemia)",
  simplify = TRUE
)

corpus_malign |>
  group_by(malignancy) |>
  tally()

corpus |>
  mutate(
    test = if_any(.cols = c(TI, AB, KW), .fns = ~ str_detect(.x, "[Aa][l][l]"))
  ) |>
  View()
# Use case-when...



# Try word cloud from "unclear" category ----------------------------------


library(wordcloud)
library(tm)

remaining_articles <- corpus_malign |>
  filter(malignancy == "Unclear")

hi <- wordcloud(remaining_articles$TI, max.words = 100)
wordcloud(remaining_articles$AB, max.words = 100)
wordcloud(remaining_articles$KW, max.words = 100)

wordcloud(corpus_malign$AB, max.words = 100)

#install.packages("wordcloud")
#install.packages("tm")

