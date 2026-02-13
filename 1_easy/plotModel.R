#!/usr/bin/env Rscript
# GSoC 2026 gfpop Easy Test: plotModel(graph) - Visualize changepoint model with time-dependent constraints
# Requires: gfpop, ggplot2

library(gfpop)
library(ggplot2)

#' Build the LOPART graph from the GSoC project description using real gfpop::Edge() calls.
#'
#' The 'rule' column is not yet supported by gfpop -- it is the proposed extension
#' for time-dependent constraints. We add it manually here to demonstrate the
#' intended API described at:
#' https://github.com/rstats-gsoc/gsoc2026/wiki/Time-dependent-constraints-in-gfpop
#'
#' @return A gfpop graph dataframe with an additional 'rule' column
create_LOPART_graph <- function() {
  # Build graph using real gfpop::Edge() calls (matches the wiki example exactly)
  LOPART.graph <- gfpop::graph(
    gfpop::Edge("normal",   "normal",   type = "null"),               # Rule 1
    gfpop::Edge("normal",   "normal",   type = "std", penalty = 5.5), # Rule 1
    gfpop::Edge("normal",   "normal",   type = "null"),               # Rule 2
    gfpop::Edge("noChange", "noChange", type = "null"),               # Rule 2
    gfpop::Edge("normal",   "noChange", type = "std"),                # Rule 2
    gfpop::Edge("normal",   "normal",   type = "std"),                # Rule 3
    gfpop::Edge("noChange", "normal",   type = "null"),               # Rule 3
    gfpop::Edge("normal",   "normal",   type = "null")                # Rule 4
  )

  # 'rule' is not yet supported by gfpop::Edge(); add it manually to demonstrate
  # the proposed time-dependent constraint API
  LOPART.graph$rule <- c(
    1, 1,      # Rule 1: unlabeled
    2, 2, 2,   # Rule 2: in positive label
    3, 3,      # Rule 3: end positive label
    4          # Rule 4: in negative label
  )

  # Semantic labels for each rule ID (used by plotModel for X-axis annotation)
  attr(LOPART.graph, "rule_labels") <- c(
    "1" = "unlabeled",
    "2" = "in positive label",
    "3" = "end positive label",
    "4" = "in negative label"
  )

  return(LOPART.graph)
}

#' Plot a changepoint model graph with time-dependent constraints
#'
#' Draws a matrix of nodes (rows = states, columns = rules) with directed edges.
#' Self-loops are drawn as curved arcs. Node labels mask edge endpoints for clarity.
#' Edge colors distinguish edge types (null = same segment, std = changepoint).
#'
#' @param graph A dataframe with state1, state2, type, penalty, and rule columns.
#'   Optionally, pass \code{rule_labels} as an attribute to show semantic labels
#'   on the X-axis (e.g., "unlabeled", "in positive label").
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

  # --- Edge type color palette ---
  # null = stay in same segment (no changepoint), std = changepoint transition
  edge_types  <- sort(unique(graph$type))
  type_colors <- c(
    "null" = "#1b9e77",   # teal -- same segment
    "std"  = "#d95f02"    # orange -- changepoint
  )
  # Fallback for any other edge types (up, down, abs, etc.)
  extra_types <- setdiff(edge_types, names(type_colors))
  if (length(extra_types) > 0) {
    extra_palette <- scales::hue_pal()(length(extra_types))
    type_colors <- c(type_colors, setNames(extra_palette, extra_types))
  }

  # --- Readable type labels for the legend ---
  type_labels <- c(
    "null" = "null (same segment)",
    "std"  = "std (changepoint)"
  )
  # Fallback for unlabeled types
  for (et in edge_types) {
    if (is.na(type_labels[et]) || is.null(type_labels[et])) {
      type_labels[et] <- et
    }
  }

  # --- Build edge labels: only show penalty when > 0 ---
  edge_label <- ifelse(
    graph$penalty > 0,
    paste0("pen=", graph$penalty),
    ""
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
  # Self-loops arc ABOVE the node: start and end are offset upward + horizontally
  # so the curve forms a clean arc above (or below for stacked loops)
  if (nrow(self_loops) > 0) {
    h_offset <- 0.18  # horizontal spread
    v_offset <- 0.18  # vertical lift so arc clears the node
    self_loops <- do.call(rbind, lapply(
      split(self_loops, interaction(self_loops$x, self_loops$y)),
      function(grp) {
        n <- nrow(grp)
        # Stack multiple self-loops: 1st above, 2nd further above, etc.
        grp$stack_idx <- seq_len(n)
        grp
      }
    ))
    self_loops$y_lift   <- v_offset + (self_loops$stack_idx - 1) * 0.22
    self_loops$x_start  <- self_loops$x - h_offset
    self_loops$y_start  <- self_loops$y + self_loops$y_lift
    self_loops$x_end    <- self_loops$x + h_offset
    self_loops$y_end    <- self_loops$y + self_loops$y_lift
    self_loops$label_y  <- self_loops$y + self_loops$y_lift + 0.15 * self_loops$stack_idx
  }

  # ================================================================
  # BUILD PLOT -- layer order: edges first, then labels on top
  # ================================================================
  p <- ggplot()

  # --- Layer 1: regular directed edges (straight arrows, colored by type) ---
  if (nrow(regular_edges) > 0) {
    # Shorten both ends so arrows don't hide under node labels
    regular_edges$y_short    <- regular_edges$y +
      0.12 * sign(regular_edges$yend - regular_edges$y)
    regular_edges$yend_short <- regular_edges$yend +
      0.12 * sign(regular_edges$y - regular_edges$yend)
    regular_edges$edge_color <- type_colors[regular_edges$type]
    p <- p + geom_segment(
      data = regular_edges,
      aes(x = x, y = y_short, xend = xend, yend = yend_short, color = type),
      arrow     = arrow(length = unit(0.2, "cm"), type = "closed"),
      linewidth = 0.8
    )
    # Show penalty label only when > 0 (white background masks line behind it)
    labeled_edges <- regular_edges[regular_edges$label != "", , drop = FALSE]
    if (nrow(labeled_edges) > 0) {
      labeled_edges$ymid <- (labeled_edges$y + labeled_edges$yend) / 2
      for (i in seq_len(nrow(labeled_edges))) {
        p <- p + geom_label(
          data = data.frame(x = labeled_edges$x[i] + 0.12,
                            y = labeled_edges$ymid[i],
                            label = labeled_edges$label[i]),
          aes(x = x, y = y, label = label),
          size          = 2.8,
          color         = labeled_edges$edge_color[i],
          fill          = "white",
          fontface      = "italic",
          linewidth     = NA,
          label.padding = unit(0.15, "lines"),
          label.r       = unit(0.1, "lines")
        )
      }
    }
  }

  # --- Layer 2: self-loop edges (curved arcs above nodes, colored by type) ---
  if (nrow(self_loops) > 0) {
    for (i in seq_len(nrow(self_loops))) {
      edge_color <- type_colors[self_loops$type[i]]
      p <- p + geom_curve(
        data = self_loops[i, ],
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
        curvature  = -0.5,
        arrow      = arrow(length = unit(0.18, "cm"), type = "closed"),
        arrow.fill = edge_color,
        color      = edge_color,
        linewidth  = 0.8
      )
    }
    # Show penalty labels for self-loops with penalty > 0 (white background)
    labeled_loops <- self_loops[self_loops$label != "", , drop = FALSE]
    if (nrow(labeled_loops) > 0) {
      for (i in seq_len(nrow(labeled_loops))) {
        edge_color <- type_colors[labeled_loops$type[i]]
        p <- p + geom_label(
          data = data.frame(x = labeled_loops$x[i],
                            y = labeled_loops$label_y[i] + 0.03,
                            label = labeled_loops$label[i]),
          aes(x = x, y = y, label = label),
          size          = 2.8,
          color         = edge_color,
          fill          = "white",
          fontface      = "italic",
          linewidth     = NA,
          label.padding = unit(0.15, "lines"),
          label.r       = unit(0.1, "lines")
        )
      }
    }
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

  # --- Color scale with legend ---
  p <- p +
    scale_color_manual(
      name   = "Transition Type",
      values = type_colors,
      labels = type_labels[names(type_colors)],
      breaks = intersect(names(type_colors), edge_types)
    )

  # --- Semantic rule labels for X-axis (if provided as attribute) ---
  rule_labs <- attr(graph, "rule_labels")
  if (is.null(rule_labs)) {
    x_labels <- paste("Rule", rules)
  } else {
    x_labels <- paste0("Rule ", rules, "\n", rule_labs[as.character(rules)])
  }

  # --- Scales, labels, theme ---
  p <- p +
    scale_x_continuous(
      breaks = seq_along(rules),
      labels = x_labels,
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
      subtitle = "Arrows show allowed transitions per rule; color distinguishes edge type"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "grey92", linetype = "dashed"),
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 10),
      axis.title.x      = element_text(face = "bold", margin = margin(t = 20)),
      axis.title.y      = element_text(face = "bold", margin = margin(r = 15)),
      axis.text         = element_text(size = 11),
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold"),
      legend.margin     = margin(t = 20),
      plot.margin       = margin(15, 15, 15, 15)
    ) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))

  return(p)
}

# Main execution: run when script is sourced or executed
run_main <- function() {
  # Build LOPART graph using real gfpop::Edge() calls + manually added 'rule' column
  LOPART.graph <- create_LOPART_graph()

  # Generate and save plot
  p <- plotModel(LOPART.graph)
  ggsave("output_plot.png", plot = p, width = 10, height = 6, dpi = 150)

  message("Plot saved to output_plot.png")
  invisible(p)
}

# Run when executed via Rscript (command line)
if (!interactive()) {
  run_main()
}
