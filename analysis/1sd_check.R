summ_all |>
  select(term, estimate, std.error, method) |>
  filter(method == "CCA") |>
  mutate(est_CCA = estimate, std_CCA = std.error) |>
  select(term, est_CCA, std_CCA) |>
  right_join(
    summ_all |>
      select(term, estimate, std.error, method) |>
      filter(!(method %in% c("CCA", "Missing\nindicator")))
  ) |>
  mutate(
    low_1std_CCA = est_CCA - std_CCA,
    upp_1std_CCA = est_CCA + std_CCA,
    within_1std = low_1std_CCA < estimate & estimate < upp_1std_CCA
  ) |>
  View()
