ggplot_grouped_forest <- function(dat,
                                  dictionary,
                                  results,
                                  event = "Relapse",
                                  form,
                                  breaks_x = c(0.5, 1, 1.5, 2, 3),
                                  lims_x = c(0.65, 4.5)) {

  # Prepare plot df
  df_plot <- prepare_forest_df(
    dat = dat,
    dictionary = dictionary,
    form = form,
    results = results
  )

  # general dodge
  dodge_w <- 0.85

  df_plot$count_all <- df_plot$count_REL + df_plot$count_NRM

  # Prepare labels and counts
  event_counts <- if (event == "Relapse") {
    "count_REL"
  } else if (event == "Non-relapse mortality") {
    "count_NRM"
  } else "count_all"

  # Prepare df for table part of forest plot
  vars_table <- c("levels_lab", "count", event_counts, "graph_x", "colour_row")
  expr_table <- parse(text = paste0(".(", paste(vars_table, collapse = ","), ")"))
  table_df <- unique(df_plot[, eval(expr_table)])
  max_x <- max(table_df$graph_x) + 2

  # Make table
  table_plot <- table_df %>%
    ggplot2::ggplot(ggplot2::aes(y = .data$levels_lab)) +
    ggplot2::geom_rect(
      xmin = -Inf,
      xmax = Inf,
      ggplot2::aes(
        fill = .data$colour_row,
        ymin = .data$graph_x - 0.5,
        ymax = .data$graph_x + 0.5
      )
    ) +
    ggplot2::geom_text(ggplot2::aes(label = .data$levels_lab), x= 0, hjust = 0, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$count), x= 1, hjust = 1, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label = !!rlang::sym(event_counts)), x = 1.3, hjust = 1, na.rm = TRUE) +
    ggplot2::annotate("text", x = 0, y = max_x, label = "Variable", hjust = 0, fontface = "bold") +
    ggplot2::annotate("text", x = 1, y = max_x, label = "n", hjust = 1, fontface = "bold") +
    ggplot2::annotate("text", x = 1.3, y = max_x, label = "# Events", hjust = 1, fontface = "bold") +
    ggplot2::coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 40)) +
    ggplot2::scale_fill_manual(values = c("gray90", "white"), guide = "none") +
    ggplot2::scale_y_discrete(limits = rev(levels(table_df$levels_lab))) +
    ggplot2::theme_void(base_size = 14)

  # Make the right hand side
  estimates_plot <- df_plot %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$levels_lab,
        y = .data$estimate,
        group = .data$method
      )
    ) +
    ggplot2::scale_y_continuous(
      name = "Hazard ratio (95% CI)",
      trans = "log",
      breaks = breaks_x
    ) +
    ggplot2::geom_rect(
      ymin = -Inf,
      ymax = Inf,
      ggplot2::aes(
        fill = .data$colour_row,
        xmin = .data$graph_x - 0.5,
        xmax = .data$graph_x + 0.5
      )
    ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = .data$`2.5 %`,
        ymax = .data$`97.5 %`,
        xmin = .data$levels_lab,
        xmax = .data$levels_lab,
        col = .data$method
      ),
      position = ggplot2::position_dodge(width = dodge_w),
      size = 0.5,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(col = .data$method, shape = .data$method),
      position = ggplot2::position_dodge(width = dodge_w),
      size = 1.5,
      na.rm = TRUE
    ) +
    ggplot2::coord_flip(ylim = lims_x, xlim = c(0, 40)) +
    ggplot2::geom_segment(x = 0, y = log(1), xend = max_x - 1, yend = log(1), linetype = "dashed") +
    ggplot2::geom_segment(x = 0, y = log(2), xend = max_x - 1, yend = log(2), linetype = "dotted") +
    ggplot2::annotate("text", x = max_x, y = 1, label = event, hjust = 0.5, fontface = "bold") +
    ggplot2::scale_x_discrete(limits = rev(levels(table_df$levels_lab))) +
    ggplot2::scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
    ggplot2::scale_shape_discrete(na.translate = FALSE) +
    ggplot2::scale_fill_manual(values = c("gray90", "white"), guide = "none") +
    ggplot2::guides(
      colour = ggplot2::guide_legend("Method", byrow = TRUE, reverse = TRUE),
      shape = ggplot2::guide_legend("Method", byrow = TRUE, reverse = TRUE)
    ) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(colour = "black", size = 12),
      axis.ticks.x = ggplot2::element_line(colour = "grey10"),
      axis.ticks.length.x = ggplot2::unit(.15, "cm"),
      axis.line.x = ggplot2::element_line(colour = "grey10"),
      axis.text.x = ggplot2::element_text(colour = "black", size = 10)
    )

  # Bring together with patchwork
  forest_plot <- table_plot +
    estimates_plot +
    patchwork::plot_layout(guides = "collect", widths = c(1, 1)) &
    ggplot2::theme(legend.position = "top", plot.margin = ggplot2::margin(0, 0, 0, 2))

  return(forest_plot)
}


prepare_forest_df <- function(dat,
                              dictionary,
                              form,
                              results) {

  # For checks
  . <- colour_row <- graph_x <- level_num <- NULL
  var_label <- var_name <- levels_lab <- term <- method  <- NULL

  # Combine list of regression results - add also NA row with reference categories
  res <- data.table::rbindlist(
    lapply(X = results, FUN = function(mod_summary) {
      reflevels_add_summary(summ = mod_summary, dat = dat, form = form)
    }),
    fill = TRUE,
    idcol = "method"
  )

  res[, method := factor(method, levels = names(results))]

  # Get predictors from RHS
  patt <- paste0("(", paste(unique(dictionary$var_name), collapse = "|"), ")")

  # Set var_name and levels according to data dict?
  # Possibly along with eventual function create_data_dictionary()
  res[, ':=' (
    var_name = regmatches(x = term, m = regexpr(pattern = patt, text = term)),
    levels = gsub(pattern = patt, replacement = "", x = term)
  )]

  res[levels == "", levels := NA_character_]
  df_merged <- merge(res, dictionary, by = c("var_name", "levels"))

  # Add rows for variable names of factors (once for each method)
  new_rows <- df_merged[!is.na(levels_lab), list(
    var_name = unique(var_name)
  ), by = list(var_label, method)]

  df_merged_expanded <- rbind(df_merged, new_rows, fill = TRUE)
  df_merged_expanded[is.na(level_num), level_num := 0]

  # Prepare the left most column of forest table, use tabs for levels
  df_merged_expanded[, levels_lab := ifelse(
    is.na(levels_lab),
    var_label,
    paste0("    ", levels_lab)
  )]

  # Sort a) alphabetic var labels, then by levels order, then method
  data.table::setorder(df_merged_expanded, var_label, level_num, method)

  # Unique() will take order they were sorted in
  df_merged_expanded[, levels_lab := factor(levels_lab, levels = unique(levels_lab))]

  # Prepare the alternating grey/white rows
  df_merged_expanded[, graph_x := as.numeric(levels_lab)]
  df_merged_expanded[order(levels_lab), colour_row := rep(
    x = c("white", "gray"), length.out = .N
  ), by = .(method)]

  return(df_merged_expanded)
}

#' Utility to add reference factor levels to add model summary
#'
#' (Not for use beyond this repository)
#'
#' @param summ Model summary (possibly custom made)
#' @param dat Original dataset used to run the model
#' @param form Formula from the run model
#' @param term_col Column referencing coefficient names in summ
#'
#' @noRd
reflevels_add_summary <- function(summ, dat, form, term_col = "term") {

  variable <- coef <- NULL

  # Get predictors
  preds <- attr(stats::terms(form), "term.labels")

  # Identify factors
  ref_levels <- dat[, lapply(.SD, function(col) levels(col)[1]), .SDcols = is.factor] %>%
    data.table::transpose(keep.names = "variable") %>%
    data.table::setnames(old = "V1", new = "coef")

  ref_levels[, "term" := paste0(variable, coef)]
  ref_levels[, setdiff(names(ref_levels), "term") := NULL]

  # Make sure term column in summary is called "term"
  new_summary <- data.table::data.table(summ)
  data.table::setnames(new_summary, old = term_col, new = "term")

  return(rbind(new_summary, ref_levels, fill = TRUE))
}
