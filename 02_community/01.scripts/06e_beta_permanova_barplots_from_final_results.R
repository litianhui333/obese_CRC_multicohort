library(tidyverse)

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
uni_dir <- file.path(base_dir, "03.results/02_beta_diversity/permanova_subgroups_univariate")
multi_dir <- file.path(base_dir, "03.results/02_beta_diversity/permanova_subgroups_multi")

uni_file <- file.path(uni_dir, "PERMANOVA_By_Continent_Results_with_Overall_remove_one_study_20260122.csv")
multi_file <- file.path(multi_dir, "PERMANOVA_Subgroup_Comparison_Summary.csv")

stopifnot(file.exists(uni_file), file.exists(multi_file))

theme_clean <- function() {
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30", size = 10),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      strip.background = element_rect(fill = "#F2F2F2", color = "#D0D0D0"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

sig_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

fill_colors <- c(
  "Disease" = "#C44E52",
  "BMI Group" = "#1B9E77",
  "BMI" = "#4C72B0",
  "Age" = "#7A7A7A",
  "Sex" = "#B8B8B8",
  "BMI X Disease" = "#8172B2",
  "BMI Group X Disease" = "#DD8452"
)

# ------------------------------------------------------------------------------
# Univariate
# ------------------------------------------------------------------------------
uni_raw <- read.csv(uni_file, check.names = FALSE, stringsAsFactors = FALSE)
if (names(uni_raw)[1] == "") {
  uni_raw <- uni_raw[, -1, drop = FALSE]
}

uni_plot_data <- uni_raw %>%
  mutate(
    Display = recode(
      Factor,
      "BMI_Group1" = "BMI Group",
      "BMI" = "BMI",
      "Disease" = "Disease",
      "Age" = "Age",
      "Sex" = "Sex"
    ),
    Display = factor(Display, levels = c("Disease", "BMI Group", "BMI", "Age", "Sex")),
    Region = factor(Region, levels = c("Overall", "Asia", "Europe")),
    R2_Percent = R2 * 100,
    Label = sig_label(P_value)
  ) %>%
  filter(!is.na(Display))

uni_clean_out <- file.path(uni_dir, "PERMANOVA_By_Continent_Results_with_Overall_remove_one_study_cleaned_20260427.csv")
write.csv(uni_plot_data, uni_clean_out, row.names = FALSE)

p_uni <- ggplot(uni_plot_data, aes(x = Display, y = R2_Percent, fill = Display)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = Label), vjust = -0.45, size = 4.5, fontface = "bold") +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Univariate PERMANOVA",
    subtitle = "Final remove-one-study results",
    y = "R2 (%)",
    caption = "* P < 0.05, ** P < 0.01, *** P < 0.001"
  ) +
  theme_clean()

ggsave(file.path(uni_dir, "PERMANOVA_R2_Plot_uni_revised_20260427.pdf"), p_uni, width = 6.6, height = 3.5)
ggsave(file.path(uni_dir, "PERMANOVA_R2_Plot_uni_revised_20260427.png"), p_uni, width = 6.6, height = 3.5, dpi = 300)

# ------------------------------------------------------------------------------
# Multivariate
# ------------------------------------------------------------------------------
multi_raw <- read.csv(multi_file, stringsAsFactors = FALSE)

multi_plot_data <- bind_rows(
  multi_raw %>%
    filter(BMI_Type == "BMI") %>%
    transmute(
      Region,
      Model = "Continuous BMI model",
      Display = "Disease",
      R2 = R2_Disease,
      P_value = P_Disease
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI") %>%
    transmute(
      Region,
      Model = "Continuous BMI model",
      Display = "BMI",
      R2 = R2_BMI,
      P_value = P_BMI
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI") %>%
    transmute(
      Region,
      Model = "Continuous BMI model",
      Display = "BMI X Disease",
      R2 = R2_Interaction,
      P_value = P_Interaction
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI_Group1") %>%
    transmute(
      Region,
      Model = "BMI-group model",
      Display = "Disease",
      R2 = R2_Disease,
      P_value = P_Disease
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI_Group1") %>%
    transmute(
      Region,
      Model = "BMI-group model",
      Display = "BMI Group",
      R2 = R2_BMI,
      P_value = P_BMI
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI_Group1") %>%
    transmute(
      Region,
      Model = "BMI-group model",
      Display = "BMI Group X Disease",
      R2 = R2_Interaction,
      P_value = P_Interaction
    )
) %>%
  mutate(
    Region = factor(Region, levels = c("Overall", "Asia", "Europe")),
    Model = factor(Model, levels = c("Continuous BMI model", "BMI-group model")),
    Display = factor(
      Display,
      levels = c("Disease", "BMI", "BMI X Disease", "BMI Group", "BMI Group X Disease")
    ),
    R2_Percent = R2 * 100,
    Label = sig_label(P_value)
  )

multi_clean_out <- file.path(multi_dir, "PERMANOVA_Subgroup_Comparison_Summary_long_cleaned_20260427.csv")
write.csv(multi_plot_data, multi_clean_out, row.names = FALSE)

multi_single_row_data <- bind_rows(
  multi_raw %>%
    group_by(Region) %>%
    summarise(
      Display = "Disease",
      R2 = mean(R2_Disease, na.rm = TRUE),
      P_value = min(P_Disease, na.rm = TRUE),
      .groups = "drop"
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI") %>%
    transmute(
      Region,
      Display = "BMI",
      R2 = R2_BMI,
      P_value = P_BMI
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI") %>%
    transmute(
      Region,
      Display = "BMI X Disease",
      R2 = R2_Interaction,
      P_value = P_Interaction
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI_Group1") %>%
    transmute(
      Region,
      Display = "BMI Group",
      R2 = R2_BMI,
      P_value = P_BMI
    ),
  multi_raw %>%
    filter(BMI_Type == "BMI_Group1") %>%
    transmute(
      Region,
      Display = "BMI Group X Disease",
      R2 = R2_Interaction,
      P_value = P_Interaction
    )
) %>%
  mutate(
    Region = factor(Region, levels = c("Overall", "Asia", "Europe")),
    Display = factor(
      Display,
      levels = c("Disease", "BMI Group", "BMI", "BMI Group X Disease", "BMI X Disease")
    ),
    R2_Percent = R2 * 100,
    Label = sig_label(P_value)
  )

multi_single_row_out <- file.path(multi_dir, "PERMANOVA_Subgroup_Comparison_Summary_single_row_cleaned_20260427.csv")
write.csv(multi_single_row_data, multi_single_row_out, row.names = FALSE)

p_multi <- ggplot(multi_single_row_data, aes(x = Display, y = R2_Percent, fill = Display)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = Label), vjust = -0.45, size = 4.2, fontface = "bold") +
  facet_wrap(~ Region, nrow = 1) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Multivariable PERMANOVA",
    subtitle = "Final adjusted results from saved summary tables",
    y = "R2 (%)",
    caption = "* P < 0.05, ** P < 0.01, *** P < 0.001"
  ) +
  theme_clean()

ggsave(file.path(multi_dir, "PERMANOVA_R2_Comparison_Final_revised_20260427.pdf"), p_multi, width = 6.8, height = 3.8)
ggsave(file.path(multi_dir, "PERMANOVA_R2_Comparison_Final_revised_20260427.png"), p_multi, width = 6.8, height = 3.8, dpi = 300)

message("Revised PERMANOVA barplots saved.")
