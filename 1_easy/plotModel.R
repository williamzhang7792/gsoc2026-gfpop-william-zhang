# plotModel.R -- GSoC 2026 Easy Test
# Visualize a changepoint model with time-dependent constraints

library(gfpop)
library(ggplot2)

# Build the LOPART graph from the GSoC wiki example using real gfpop calls.
# The 'rule' column doesn't exist in gfpop yet -- that's the proposed extension.
# We add it manually here to show what the API would look like.
create_LOPART_graph <- function() {
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

  # Add rule column manually (not yet in gfpop)
  LOPART.graph$rule <- c(
    1, 1,      # Rule 1: unlabeled
    2, 2, 2,   # Rule 2: in positive label
    3, 3,      # Rule 3: end positive label
    4          # Rule 4: in negative label
  )

  # Semantic labels for X-axis
  attr(LOPART.graph, "rule_labels") <- c(
    "1" = "unlabeled",
    "2" = "in positive label",
    "3" = "end positive label",
    "4" = "in negative label"
  )

  return(LOPART.graph)
}

# Plot a gfpop graph with time-dependent constraints.
# Rows = states, columns = rules, edges = allowed transitions.
plotModel <- function(graph) {
  if (!"rule" %in% names(graph))
    stop("Graph must contain a 'rule' column.")

  states <- sort(unique(c(graph$state1, graph$state2)))
  rules  <- sort(unique(graph$rule))

  state_to_y <- setNames(seq_along(states), states)
  rule_to_x  <- setNames(seq_along(rules), rules)

  edge_types  <- sort(unique(graph$type))
  type_colors <- c("null" = "#1b9e77", "std" = "#d95f02")
  type_labels <- c("null" = "null (same segment)", "std" = "std (changepoint)")

  # Only show penalty when > 0
  edge_label <- ifelse(
    graph$penalty > 0,
    paste0("pen=", graph$penalty),
    ""
  )

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

  node_data <- expand.grid(rule = rules, state = states, stringsAsFactors = FALSE)
  node_data$x <- rule_to_x[as.character(node_data$rule)]
  node_data$y <- state_to_y[as.character(node_data$state)]

  # Self-loops: arc above the node with horizontal + vertical offset
  if (nrow(self_loops) > 0) {
    h_offset <- 0.18
    v_offset <- 0.18
    self_loops <- do.call(rbind, lapply(
      split(self_loops, interaction(self_loops$x, self_loops$y)),
      function(grp) {
        grp$stack_idx <- seq_len(nrow(grp))
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

  # Build plot: edges first, then nodes on top
  p <- ggplot()

  # Straight arrows (shorten ends so arrowheads don't hide under labels)
  if (nrow(regular_edges) > 0) {
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

  # Self-loop arcs
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

  # Nodes on top (white fill masks the line endings)
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

  p <- p +
    scale_color_manual(
      name   = "Transition Type",
      values = type_colors,
      labels = type_labels[names(type_colors)],
      breaks = intersect(names(type_colors), edge_types)
    )

  rule_labs <- attr(graph, "rule_labels")
  if (is.null(rule_labs)) {
    x_labels <- paste("Rule", rules)
  } else {
    x_labels <- paste0("Rule ", rules, "\n", rule_labs[as.character(rules)])
  }

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

# Run
LOPART.graph <- create_LOPART_graph()
p <- plotModel(LOPART.graph)
ggsave("output_plot.png", plot = p, width = 10, height = 6, dpi = 150)
message("Plot saved to output_plot.png")
