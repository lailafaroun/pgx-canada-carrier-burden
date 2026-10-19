library(tidyverse)

dir.create("figures", showWarnings = FALSE)

# Shared project palette: yellow -> orange -> red
pgx_colours <- c(
  "#FFF7BC",
  "#FEC44F",
  "#FD8D3C",
  "#F03B20",
  "#BD0026"
)

df <- read_csv("data/carrier_by_group.csv", show_col_types = FALSE)

required_columns <- c(
  "group", "gene", "carrier", "expresser", "low_an", "few_alt"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing expected columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

group_order <- c(
  "African/African American",
  "Admixed American",
  "Ashkenazi Jewish",
  "East Asian",
  "European (Finnish)",
  "European (non-Finnish)",
  "Middle Eastern",
  "Remaining",
  "South Asian",
  "Amish"
)

gene_order <- c(
  "ABCG2",
  "CYP3A5",
  "CYP4F2",
  "NUDT15",
  "CYP3A5 expresser",
  "VKORC1"
)

# Main variant carrier estimates
carrier_data <- df %>%
  transmute(
    group,
    gene,
    estimate = carrier,
    flagged = coalesce(low_an, FALSE) |
      coalesce(few_alt, FALSE)
  )

# CYP3A5 expresser estimates are a separate measure
expresser_data <- df %>%
  filter(gene == "CYP3A5", !is.na(expresser)) %>%
  transmute(
    group,
    gene = "CYP3A5 expresser",
    estimate = expresser,
    flagged = coalesce(low_an, FALSE) |
      coalesce(few_alt, FALSE)
  )

plot_data <- bind_rows(carrier_data, expresser_data) %>%
  mutate(
    percent = estimate * 100,
    label = paste0(
      round(percent),
      if_else(flagged, "*", "")
    ),
    group = factor(group, levels = rev(group_order)),
    gene = factor(gene, levels = gene_order),
    label_colour = if_else(percent >= 70, "white", "#222222")
  )

if (any(is.na(plot_data$group))) {
  stop("Unexpected ancestry-group names found in the input CSV.")
}

p <- ggplot(plot_data, aes(x = gene, y = group, fill = percent)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = label, colour = label_colour),
    size = 3.2,
    show.legend = FALSE
  ) +
  scale_fill_gradientn(
    colours = pgx_colours,
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    name = "Probability (%)",
    oob = scales::squish
  ) +
  scale_colour_identity() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "Population-level pharmacogenomic variant burden",
    subtitle = "Estimated carrier probabilities across gnomAD ancestry groups",
    x = NULL,
    y = NULL,
    caption = paste(
      "* Low data support: allele number (AN) < 2,500",
      "or alternate allele count (AC) < 20.",
      "CYP3A5 expresser probability is shown separately",
      "from CYP3A5*3 carrier probability.",
      "Combined exome and genome counts are used where available."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 35, hjust = 1, colour = "#222222"
    ),
    axis.text.y = element_text(colour = "#222222"),
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      colour = "#444444"
    )
  )

ggsave(
  "figures/fig1_carrier_heatmap.png",
  plot = p,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/fig1_carrier_heatmap.pdf",
  plot = p,
  width = 11,
  height = 7,
  bg = "white"
)

cat("Saved Figure 1 PNG and PDF in figures/\\n")
