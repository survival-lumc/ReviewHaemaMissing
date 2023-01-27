ggplot_grouped_forest <- function(dat,
                                  dictionary,
                                  results,
                                  event = "Relapse",
                                  form,
                                  breaks_x = c(0.5, 1, 1.5, 2, 3),
                                  lims_x = c(0.65, 4.5)) {

  # Prepare plot df
  df_plot <- CauseSpecCovarMI:::prepare_forest_df(
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
