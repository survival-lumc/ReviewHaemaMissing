# Load libraries
library(readxl)
library(tidyverse)

# Read-in sheet
corpus_raw <- readxl::read_xls(path = "../2022-07-06_export.xls", skip = 1L)

# Remove useless columns
corpus <- corpus_raw |>
  mutate(
    AU = str_replace_all(AU, pattern = "(\n\n)", replacement = ","),
    KW = str_replace_all(KW, pattern = "(\n\n)", replacement = ",")
  ) |>
  select(-all_of(c("VN", "FTURL", "XL"))) |>

  # Clean journal names - cross-check with IO
  mutate(JA = ifelse(is.na(JA), "Transplant. Cell. Ther", JA)) |>

  # Reorder
  relocate(c("AB", "KW"), .after = c("TI", "AU"))

# What journals?
corpus |>
  ggplot(aes(reorder(JA, JA, length))) + #fct_infreq also possible
  geom_bar(fill = viridisLite::viridis(1)) +
  stat_count(
    geom = "text", colour = "white", size = 3.5,
    aes(label = ..count..),position=position_stack(vjust=0.95)
  ) +
  coord_flip() +
  theme_minimal()

# EBMT authors
corpus |> filter(str_detect(AU, "Wreede")) # |> View() # Remove eventually?
corpus |> filter(str_detect(AU, "Schetelig"))

# Checking author clustering
corpus |>
  mutate(n_total_articles = n()) |>
  separate_rows(AU, sep = ",") |>
  group_by(AU, n_total_articles) |>
  tally() |>
  ungroup() |>
  mutate(prop = paste0(round(100 * n / n_total_articles, 2), "%")) |>
  filter(n >= 5) |>
  arrange(desc(n))# |>
  #View()

# Check clustering of CIBMTR/EBMT -> check in abstract too
filter(.data = corpus, if_any(c("TI", "AB"), ~ str_detect(.x, pattern = "EBMT")))
filter(.data = corpus, if_any(c("TI", "AB"), ~ str_detect(.x, pattern = "CIBMTR")))


# Disease regexes ---------------------------------------------------------

regex_malign <- list(
  "AML" = "AML|([Aa]cute [Mm]yeloid [Ll]euk[ae]*)|([Aa]cute [Mm]yelogenous [Ll]euk[ae]*)",
  "NHL" = "[Nn]on-Hodgkin*",
  "MDS" = "([Mm]yelodysplastic [Ss]yndrome)",
  "CLL" = "([Cc]hronic [Ll]ymphocytic [Ll]euk[ae]*)",
  "CML" = "([Cc]hronic [Mm]yeloid [Ll]euk[ae]*)|([Cc]hronic [Mm]yelogenous [Ll]euk[ae]*)",
  "ALL" = "([Aa]cute [Ll]ymphocytic [Ll]euk[ae]*)|([Aa]cute [Ll]ymphoblastic [Ll]euk[ae]*)",
  "MM" = "([Mm]ultiple [Mm]yeloma*)|[Mm]yeloma*",
  "AL" = "([Aa]cute leuk[ae]*)", # No clue if AML or AL
  "APL" = "([Aa]cute [Pp]romyelocytic leuk[ae]*)",
  "leukem_other_raw" = "leuk[ae]*"
)

# The master regex - if needed
regex_master <- str_flatten(regex_malign, "|")


# Overview ----------------------------------------------------------------


# Add "all_fields" to search
corpus <- corpus |>
  mutate(all_fields = paste(TI, AB, KW, sep = " ; "))

# Time to classify papers
malign_df_raw <- data.frame(
  sapply(regex_malign, function(reg) as.integer(str_detect(corpus$all_fields, pattern = reg)))
)

# Keep leukem_other as true only if not used at same time as AML, CLL, CML, ALL
malign_df <- malign_df_raw |>
  mutate(
    leukem_other = ifelse(if_any(c("AML", "CLL", "CML", "ALL", "AL", "APL"), ~ .x == 1), 0L, leukem_other_raw)
  ) |>
  select(-leukem_other_raw) |>
  mutate(n_malign = rowSums(across(everything())))

# Four of them have none - these are ones using "mm" as unit or multiple / things of the sort
table(malign_df$n_malign)

# Noice
UpSetR::upset(data = malign_df |> select(-n_malign), nsets = 10)

# Checking the Labopin cluster --------------------------------------------

malign_df |>
  slice(which(str_detect(corpus$AU, "Labopin"))) |>
  select(-n_malign) |>
  UpSetR::upset(nsets = 10)
#
#corpus |> filter() |>  View()

# Check out the myeloma
corpus |>
  slice(which(with(malign_df, MM == 1 & AML == 1))) |>
  slice(which(str_detect(AU, "Labopin"))) |>
  View()

# Have a look at those "other leukemias" ----------------------------------

# Specifically, those used without any other combo

corpus |>
  slice(which(with(malign_df, leukem_other == 1 & n_malign == 1))) |>
  View()

# What are we finding? Use master regex
corpus |>
  slice(which(with(malign_df, leukem_other == 1 & n_malign == 1))) |>
  pull(all_fields) |>
  wordcloud::wordcloud(max.words = 100)


# .. some acute leukemias (AL), some specific ones like Burkitt,
# some also myelofibrosis (since it can transform to Leukemia!)



# Try and remove ALL ------------------------------------------------------


# Note how do we deal with acute leukemia (myeloid vs lympocytics)

corpus |>
  slice(
    -which(with(malign_df, (ALL == 1 | NHL == 1 | CLL == 1) & n_malign == 1))
  ) |>
  mutate(n_total_articles = n()) |>
  separate_rows(AU, sep = ",") |>
  group_by(AU, n_total_articles) |>
  tally() |>
  ungroup() |>
  mutate(prop = paste0(round(100 * n / n_total_articles, 2), "%")) |>
  filter(n >= 5) |>
  arrange(desc(n))

# If we remove Labopin - no cluster with statisticians?


# Any with imp? -----------------------------------------------------------



corpus |>
  filter(str_detect(all_fields, "missing|incomplete|imput*")) |>
  View()

corpus |>
  filter(str_detect(all_fields, "missing|incomplete|imput*")) |>
  pull(all_fields) |>
  str_extract_all(pattern = "missing|incomplete|imput*")



# Check the lancet ones ---------------------------------------------------


corpus |>
  filter(JA == "Lancet Haematol.")


# need aHR in search?
