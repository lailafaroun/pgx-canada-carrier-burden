library(tidyverse)

dir.create("figures", showWarnings = FALSE)

# Shared palette used across all project figures
pgx_colours <- c(
  "#FFF7BC",
  "#FEC44F",
  "#FD8D3C",
  "#F03B20",
  "#BD0026"
)

df <- read_csv(
  "data/evidence_counts.csv",
  show_col_types = FALSE
)

required_columns <- c(
  "group",
  "unique_pmid_group_pairs",
  "group_type"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing expected columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

plot_data <- df %>%
  mutate(
    group = recode(
      group,
      "African American/Afro-Caribbean" =
        "African American /\nAfro-Caribbean",
      "Central/South Asian" = "Central /\nSouth Asian",
      "Non-population / in vitro" =
        "Non-population /\nin vitro"
    ),
    category = case_when(
      group == "Unknown" ~ "Unknown",
      group == "Multiple groups" ~ "Multiple groups",
      group_type == "Non-population" ~ "Non-population",
      TRUE ~ "Population label"
    ),
    group = fct_reorder(
      group,
      unique_pmid_group_pairs
    )
  )

# Colour categories while retaining the common project palette
category_colours <- c(
  "Population label" = pgx_colours[2],
  "Multiple groups" = pgx_colours[3],
  "Unknown" = pgx_colours[5],
  "Non-population" = pgx_colours[1]
)

p <- ggplot(
  plot_data,
  aes(
    x = unique_pmid_group_pairs,
    y = group,
    fill = category
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = unique_pmid_group_pairs),
    hjust = -0.25,
    size = 3.5,
    colour = "#222222"
  ) +
  scale_fill_manual(
    values = category_colours,
    breaks = c(
      "Population label",
      "Multiple groups",
      "Unknown",
      "Non-population"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.12))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Population representation in pharmacogenomic evidence",
    subtitle = paste(
      "Unique PMID × population-label pairs",
      "in ClinPGx drug annotations for five genes"
    ),
    x = "Number of unique PMID × population-label pairs",
    y = NULL,
    caption = paste(
      strwrap(
        paste(
          "Counts use unique PMID × population-label pairs;",
          "a PMID may contribute to multiple population labels.",
          "Unknown and Multiple groups are retained as separate categories.",
          "These are curated ClinPGx drug annotations, not a systematic review",
          "or a complete count of all published pharmacogenomic research."
        ),
        width = 100
      ),
      collapse = "\n"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "#222222"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      colour = "#444444",
      lineheight = 1.2,
      margin = margin(t = 10)
    ),
    plot.caption.position = "plot",
    plot.margin = margin(10, 25, 20, 10)
  )

ggsave(
  "figures/fig3_evidence_diversity.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/fig3_evidence_diversity.pdf",
  plot = p,
  width = 10,
  height = 8,
  bg = "white"
)

cat("Saved Figure 3 PNG and PDF in figures/\n")

print(
  plot_data %>%
    arrange(desc(unique_pmid_group_pairs)) %>%
    select(group, unique_pmid_group_pairs, category)
)
