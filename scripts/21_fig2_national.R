library(tidyverse)

dir.create("figures", showWarnings = FALSE)

# Shared palette used across all project figures
pgx_colours <- c(
  "#FFF7BC",  # pale yellow
  "#FEC44F",  # yellow
  "#FD8D3C",  # orange
  "#F03B20",  # red-orange
  "#BD0026"   # dark red
)

quantity_order <- c(
  "ABCG2",
  "CYP3A5",
  "CYP3A5_expressers",
  "CYP4F2",
  "NUDT15",
  "VKORC1"
)

display_labels <- c(
  ABCG2 = "ABCG2",
  CYP3A5 = "CYP3A5",
  CYP3A5_expressers = "CYP3A5 expresser probability",
  CYP4F2 = "CYP4F2",
  NUDT15 = "NUDT15",
  VKORC1 = "VKORC1"
)

# Read source data
scenarios <- read_csv(
  "data/all_scenario_estimates.csv",
  show_col_types = FALSE
)

sensitivity <- read_csv(
  "data/sensitivity_summary.csv",
  show_col_types = FALSE
)

# Primary national estimates and 95% Wilson intervals from S0
primary <- scenarios %>%
  filter(
    scenario == "S0",
    quantity %in% quantity_order
  ) %>%
  select(
    quantity,
    estimate = national,
    low,
    high
  )

# Sensitivity range excludes S1 because it is not nationally representative
sens <- sensitivity %>%
  filter(quantity %in% quantity_order) %>%
  select(
    quantity,
    sensitivity_low = min_excluding_S1,
    sensitivity_high = max_excluding_S1
  )

plot_data <- primary %>%
  inner_join(sens, by = "quantity") %>%
  mutate(
    quantity = factor(quantity, levels = rev(quantity_order)),
    label = recode(as.character(quantity), !!!display_labels),
    estimate_pct = estimate * 100,
    low_pct = low * 100,
    high_pct = high * 100,
    sensitivity_low_pct = sensitivity_low * 100,
    sensitivity_high_pct = sensitivity_high * 100
  )

if (nrow(plot_data) != length(quantity_order)) {
  stop(
    "Expected six matched estimates; found ",
    nrow(plot_data),
    ". Check the input CSV files."
  )
}

# Build plot
p <- ggplot(plot_data, aes(y = quantity)) +

  # Sensitivity range: orange
  geom_segment(
    aes(
      x = sensitivity_low_pct,
      xend = sensitivity_high_pct,
      yend = quantity,
      colour = "Sensitivity range"
    ),
    linewidth = 5,
    lineend = "round"
  ) +

  # 95% Wilson interval: red
  geom_segment(
    aes(
      x = low_pct,
      xend = high_pct,
      yend = quantity,
      colour = "95% Wilson interval"
    ),
    linewidth = 1.2
  ) +

  # Primary S0 estimate: yellow point
  geom_point(
    aes(
      x = estimate_pct,
      fill = "Primary estimate (S0)"
    ),
    shape = 21,
    size = 3.8,
    stroke = 0.8,
    colour = "#222222"
  ) +

  scale_colour_manual(
    name = NULL,
    values = c(
      "Sensitivity range" = pgx_colours[3],
      "95% Wilson interval" = pgx_colours[5]
    ),
    breaks = c(
      "Primary estimate (S0)",
      "95% Wilson interval",
      "Sensitivity range"
    )
  ) +

  scale_fill_manual(
    name = NULL,
    values = c("Primary estimate (S0)" = pgx_colours[2])
  ) +

  scale_y_discrete(
    labels = display_labels[levels(plot_data$quantity)]
  ) +

  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +

  labs(
    title = "Estimated pharmacogenomic variant burden in Canada",
    subtitle = paste(
      "Primary national estimates with 95% Wilson intervals",
      "and sensitivity ranges"
    ),
    x = "Estimated probability",
    y = NULL,
    caption = paste(
      strwrap(
        paste(
          "Sensitivity ranges include S0 and S2-S5; S1 is excluded because it does not represent the national population.",
          "CYP3A5 expresser probability is distinct from CYP3A5*3 carrier probability.",
          "These are modeled probabilities, not the proportion of people who will necessarily experience an adverse drug reaction or require a medication change."
        ),
        width = 100
      ),
      collapse = "\n"
    )
  ) +

  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(shape = 21, colour = "#222222", size = 3.8)
    ),
    colour = guide_legend(
      order = 2,
      override.aes = list(
        linewidth = c(1.2, 5)
      )
    )
  ) +

  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text = element_text(colour = "#222222"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 9),
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      colour = "#444444",
      lineheight = 1.2,
      margin = margin(t = 10)
    ),
    plot.caption.position = "plot",
    plot.margin = margin(10, 15, 25, 15)
  )

ggsave(
  "figures/fig2_national_estimates.png",
  plot = p,
  width = 10,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/fig2_national_estimates.pdf",
  plot = p,
  width = 10,
  height = 8.5,
  bg = "white"
)

cat("Saved Figure 2 PNG and PDF in figures/\n")

print(
  plot_data %>%
    transmute(
      quantity = label,
      estimate_pct = round(estimate_pct, 2),
      CI_low_pct = round(low_pct, 2),
      CI_high_pct = round(high_pct, 2),
      sensitivity_low_pct = round(sensitivity_low_pct, 2),
      sensitivity_high_pct = round(sensitivity_high_pct, 2)
    )
)
