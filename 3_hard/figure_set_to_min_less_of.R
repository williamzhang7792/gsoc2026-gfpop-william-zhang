# figure_set_to_min_less_of.R
# Generates a plot showing the geometric action of set_to_min_less_of:
# the isotonic constraint operator traces the input cost downhill,
# then flatlines at the minimum.

library(ggplot2)

mu <- seq(-1, 10, length.out = 500)

# Input: a simple quadratic cost C(mu) = (mu - 4)^2 + 1
cost_input <- (mu - 4)^2 + 1

# Output: Q(mu) = min_{mu' <= mu} C(mu')
# Follows the parabola down to argmin=4, then flat at cost=1.
cost_output <- ifelse(mu <= 4, cost_input, 1)

df <- data.frame(
  mu   = rep(mu, 2),
  cost = c(cost_input, cost_output),
  curve = rep(c("C(\u03bc)  [input]", "Q(\u03bc) = min C(\u03bc\u2019)  [output]"), each = length(mu))
)
df$curve <- factor(df$curve,
  levels = c("C(\u03bc)  [input]", "Q(\u03bc) = min C(\u03bc\u2019)  [output]"))

p <- ggplot(df, aes(x = mu, y = cost, color = curve, linetype = curve)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  coord_cartesian(xlim = c(-1, 10), ylim = c(0, 14)) +
  geom_vline(xintercept = 4, linetype = "dotted", color = "grey50", linewidth = 0.5) +
  annotate("point", x = 4, y = 1, color = "#D55E00", size = 3.5) +
  annotate("segment", x = 4, xend = 4, y = 1.25, yend = 3.0,
           color = "grey50", linewidth = 0.4, linetype = "solid") +
  annotate("label", x = 4, y = 3.3, label = expression(mu^"*" == 4),
           color = "grey30", hjust = 0.5, size = 3.8,
           fill = "white", linewidth = 0, label.padding = unit(0.2, "lines")) +
  annotate("segment", x = 6.8, xend = 7.5, y = 5.0, yend = 1.3,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           color = "#D55E00", linewidth = 0.5) +
  annotate("text", x = 4.5, y = 6.0,
           label = "flat: cost fixed at min\n(non-decreasing constraint)",
           color = "#D55E00", hjust = 0, size = 3.3) +
  labs(
    title = expression("Isotonic Constraint Operator: " * italic("set_to_min_less_of") * "()"),
    subtitle = expression("Enforces " * mu[k] >= mu[k-1] * " by sweeping a running minimum"),
    x = expression(mu),
    y = "Cost",
    color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = c(0.22, 0.85),
    legend.background = element_rect(fill = "white", color = "grey85", linewidth = 0.3),
    legend.margin = margin(4, 8, 4, 8),
    panel.grid.minor = element_blank()
  )

ggsave("set_to_min_less_of.png", p, width = 6, height = 4, dpi = 200)
cat("Saved set_to_min_less_of.png\n")
