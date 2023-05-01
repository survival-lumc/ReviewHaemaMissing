# Read-in raw corpus from Zotero (after importing .ris file from OVID)
corpus_raw <- read_csv("data-raw/literature-database-raw.csv")

# Re-format
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
    ~ c(
      "authors",
      "title_article",
      "abstract",
      "journal_abbrev",
      "journal_name",
      "ISSN",
      "tags",
      "DOI"
    )
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

# Simplify to first author et al.
corpus_to_export <- corpus |>
  separate_rows(authors, sep = "; ", ) |>
  group_by(title_article) |>
  slice(1L) |>
  mutate(authors = paste0(authors, " et al.")) |>
  select(authors, title_article, DOI, journal_name) |>
  arrange(authors, title_article)

# Write to extraction sheet file:
#writexl::write_xlsx(corpus_to_export, path = "data-raw//extraction-sheet.xlsx")
