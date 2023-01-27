# Decide whether to use ovid export or zotero export; here we use zotero one
library(tidyverse)

# ... need to try this now with zotero export
corpus_raw <- read_csv("data-raw/literature-database-raw.csv")

corpus <- corpus_raw |>
  select(
    Author,
    Title,
    `Abstract Note`,
    `Journal Abbreviation`,
    `Publication Title`,
    ISSN,
    `Manual Tags`,
    DOI
  ) |>
  rename_with(
    ~ c("authors", "title_article", "abstract",
        "journal_abbrev", "journal_name", "ISSN", "tags", "DOI")
  ) |>
  mutate(
    journal_abbrev = ifelse(
      journal_abbrev == "Transplant Cell Ther",
      "Transplant. Cell. Ther",
      journal_abbrev
    ),
    journal_name = ifelse(
      journal_name  == "Transplantation and cellular therapy",
      "Transplantation and Cellular Therapy",
      journal_name
    )
  )

# Remember to exlucde based on supervisors??
table(corpus$journal_abbrev) #|>  length()
table(corpus$journal_name)  #|>  length()

# Journals
corpus |>
  ggplot(aes(reorder(journal_abbrev, journal_abbrev, length))) + #fct_infreq also possible
  geom_bar(fill = viridisLite::viridis(1)) +
  stat_count(
    geom = "text", colour = "white", size = 3.5,
    aes(label = ..count..),position=position_stack(vjust=0.95)
  ) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Journal")

corpus |>
  mutate(n_total_articles = n()) |>
  separate_rows(authors, sep = "; ") |>
  group_by(authors, n_total_articles) |>
  tally() |>
  ungroup() |>
  mutate(prop = paste0(round(100 * n / n_total_articles, 2), "%")) |>
  filter(n >= 5) |>
  arrange(desc(n))

# Prep export
corpus_to_export <- corpus |>
  separate_rows(authors, sep = "; ", ) |>
  group_by(title_article) |>
  slice(1L) |>
  mutate(authors = paste0(authors, " et al.")) |>
  select(authors, title_article, DOI, journal_name) |>
  arrange(authors, title_article)

writexl::write_xlsx(corpus_to_export, path = "../review-extraction/extraction-sheet.xlsx")
