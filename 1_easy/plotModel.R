#!/usr/bin/env Rscript
# GSoC 2026 gfpop Easy Test: plotModel(graph) - Visualize changepoint model with time-dependent constraints
# Requires: ggplot2

library(ggplot2)

#' Create mock graph dataframe mimicking gfpop::graph output with added 'rule' column
#' @return A dataframe with columns: state1, state2, type, parameter, penalty, K, a, min, max, rule
create_mock_graph <- function() {
  # gfpop graph columns: state1, state2, type, parameter, penalty, K, a, min, max
  # We add 'rule' for time-dependent constraints (not yet in gfpop)
  graph <- data.frame(
    state1 = c(
      "normal", "normal",           # Rule 1
      "normal", "noChange", "normal", # Rule 2
      "normal", "noChange",          # Rule 3
      "normal"                       # Rule 4
    ),
    state2 = c(
      "normal", "normal",           # Rule 1
      "normal", "noChange", "noChange", # Rule 2
      "normal", "normal",            # Rule 3
      "normal"                       # Rule 4
    ),
    type = c(
      "null", "std",                # Rule 1
      "null", "null", "std",        # Rule 2
      "std", "null",                # Rule 3
      "null"                         # Rule 4
    ),
    parameter = 0,
    penalty = c(
      0, 5.5,    # Rule 1: null (0), std (5.5)
      0, 0, 0,   # Rule 2
      0, 0,      # Rule 3
      0          # Rule 4
    ),
    K = Inf,
    a = 0,
    min = -Inf,
    max = Inf,
    rule = c(
      1, 1,      # Rule 1
      2, 2, 2,   # Rule 2
      3, 3,      # Rule 3
      4          # Rule 4
    ),
    stringsAsFactors = FALSE
  )
  class(graph) <- c("graph", class(graph))
  return(graph)
}

#' Plot a changepoint model graph with time-dependent constraints
#'
#' Draws a matrix of nodes (rows = states, columns = rules) with directed edges.
#' Self-loops are drawn as curved arcs. Node labels mask edge endpoints for clarity.
#'
#' @param graph A dataframe with state1, state2, type, penalty, and rule columns
#' @return A ggplot object
plotModel <- function(graph) {
  if (!"rule" %in% names(graph)) {
    stop("Graph must contain a 'rule' column for time-dependent constraints.")
  }

  # --- Coordinate mapping: states -> Y-axis, rules -> X-axis ---
  states <- sort(unique(c(graph$state1, graph$state2)))
  rules  <- sort(unique(graph$rule))

  state_to_y <- setNames(seq_along(states), states)
  rule_to_x  <- setNames(seq_along(rules), rules)

  # --- Build edge labels (e.g., "std\npen=5.5") ---
  edge_label <- ifelse(
    graph$penalty > 0,
    paste0(graph$type, "\npen=", graph$penalty),
    graph$type
  )

  # --- Prepare edge data with mapped coordinates ---
  edges <- data.frame(
    x        = rule_to_x[as.character(graph$rule)],
    y        = state_to_y[as.character(graph$state1)],
    xend     = rule_to_x[as.character(graph$rule)],
    yend     = state_to_y[as.character(graph$state2)],
    type     = graph$type,
    penalty  = graph$penalty,
    label    = edge_label,
    is_self  = graph$state1 == graph$state2,
    stringsAsFactors = FALSE
  )

  regular_edges <- edges[!edges$is_self, , drop = FALSE]
  self_loops    <- edges[edges$is_self, , drop = FALSE]

  # --- Node data: full state x rule matrix ---
  node_data <- expand.grid(rule = rules, state = states, stringsAsFactors = FALSE)
  node_data$x <- rule_to_x[as.character(node_data$rule)]
  node_data$y <- state_to_y[as.character(node_data$state)]

  # --- Prepare self-loop geometry ---
  # Offset endpoints horizontally and alternate curvature for stacked loops
  if (nrow(self_loops) > 0) {
    loop_offset <- 0.35
    # Assign alternating direction per (x, y) group so stacked loops fan out
    self_loops <- do.call(rbind, lapply(
      split(self_loops, interaction(self_loops$x, self_loops$y)),
      function(grp) {
        n <- nrow(grp)
        grp$curve_dir <- rep(c(1, -1), length.out = n)
        grp
      }
    ))
    self_loops$x_start <- self_loops$x - loop_offset
    self_loops$y_start <- self_loops$y
    self_loops$x_end   <- self_loops$x + loop_offset
    self_loops$y_end   <- self_loops$y
    # Label position: above or below depending on curve direction
    self_loops$label_y_nudge <- 0.35 * self_loops$curve_dir
  }

  # ================================================================
  # BUILD PLOT -- layer order: edges first, then labels on top
  # ================================================================
  p <- ggplot()

  # --- Layer 1: regular directed edges (straight arrows) ---
  if (nrow(regular_edges) > 0) {
    # Shorten segments so arrowheads stop just outside the target node label
    regular_edges$yend_short <- regular_edges$yend +
      0.15 * sign(regular_edges$y - regular_edges$yend)
    p <- p + geom_segment(
      data = regular_edges,
      aes(x = x, y = y, xend = xend, yend = yend_short),
      arrow     = arrow(length = unit(0.25, "cm"), type = "closed"),
      arrow.fill = "grey30",
      color      = "grey30",
      linewidth  = 0.8
    )
    # Edge type labels on regular edges (midpoint)
    regular_edges$xmid <- (regular_edges$x + regular_edges$xend) / 2
    regular_edges$ymid <- (regular_edges$y + regular_edges$yend) / 2
    p <- p + geom_label(
      data = regular_edges,
      aes(x = xmid + 0.15, y = ymid, label = label),
      size       = 2.5,
      fill       = "white",
      color      = "grey50",
      label.r    = unit(0.15, "lines"),
      linewidth  = 0.2
    )
  }

  # --- Layer 2: self-loop edges (curved arrows) ---
  if (nrow(self_loops) > 0) {
    for (i in seq_len(nrow(self_loops))) {
      p <- p + geom_curve(
        data = self_loops[i, ],
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
        curvature  = 0.6 * self_loops$curve_dir[i],
        arrow      = arrow(length = unit(0.25, "cm"), type = "closed"),
        arrow.fill = "grey30",
        color      = "grey30",
        linewidth  = 0.8
      )
    }
    # Edge type labels on self-loops (nudged above/below)
    p <- p + geom_label(
      data = self_loops,
      aes(x = x, y = y + label_y_nudge, label = label),
      size       = 2.5,
      fill       = "white",
      color      = "grey50",
      label.r    = unit(0.15, "lines"),
      linewidth  = 0.2
    )
  }

  # --- Layer 3 (top): node labels, masking edge endpoints ---
  p <- p + geom_label(
    data = node_data,
    aes(x = x, y = y, label = state),
    size          = 4,
    fontface      = "bold",
    fill          = "white",
    color         = "steelblue",
    label.r       = unit(0.3, "lines"),
    linewidth     = 1.2,
    label.padding = unit(0.3, "lines")
  )

  # --- Scales, labels, theme ---
  p <- p +
    scale_x_continuous(
      breaks = seq_along(rules),
      labels = paste("Rule", rules),
      expand = expansion(mult = 0.2)
    ) +
    scale_y_continuous(
      breaks = seq_along(states),
      labels = states,
      expand = expansion(mult = 0.25)
    ) +
    labs(
      x        = "Rule (time interval)",
      y        = "State",
      title    = "Changepoint Model: Time-Dependent Constraints",
      subtitle = "Arrows indicate allowed transitions (Box shows: cost type & penalty)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linetype = "dashed"),
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 10),
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(size = 11),
      plot.margin       = margin(15, 15, 15, 15)
    )

  return(p)
}

# Main execution: run when script is sourced or executed
run_main <- function() {
  # Create mock graph data
  mock_graph <- create_mock_graph()

  # Generate and save plot
  p <- plotModel(mock_graph)
  ggsave("output_plot.png", plot = p, width = 10, height = 6, dpi = 150)

  message("Plot saved to output_plot.png")
  invisible(p)
}

# Run when executed via Rscript (command line)
if (!interactive()) {
  run_main()
}
