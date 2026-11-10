rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

base_dir <- "G:/My Drive/proj_ob_AA_CRC"
results_dir <- file.path(
  base_dir,
  "03.results/06_network_analysis/DGCA_Results/Union_DAonly_ObeseDA_Participating_CTRL_CRC_Rewiring_p01"
)
heatmap_dir <- file.path(
  base_dir,
  "03.results/04_differential_abundance/heatmap_final/BMI_Stratified_Disease_Heatmaps"
)

edge_file <- file.path(
  results_dir,
  "Union_DAonly_ObeseDA_Participating_Rewiring_Edges_Excluding_NwSameTrajectory_p01.csv"
)
out_dir <- file.path(results_dir, "Cytoscape_Excluding_NwSameTrajectory_FullGainLoss")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_files <- c(
  Normal = file.path(heatmap_dir, "Feature_List_Normal_Refined.csv"),
  Ow = file.path(heatmap_dir, "Feature_List_Ow_Refined.csv"),
  Ob = file.path(heatmap_dir, "Feature_List_Ob_Refined.csv")
)

stopifnot(file.exists(edge_file), all(file.exists(feature_files)))

format_sgb_names <- function(names_vec) {
  sapply(names_vec, function(x) {
    clean_name <- gsub(".*s__", "s__", x)
    parts <- strsplit(clean_name, "\\|")[[1]]
    sp_part <- gsub("^s__", "", parts[1])
    sgb_id <- if (length(parts) > 1) gsub("^t__", "", parts[2]) else ""

    if (grepl("_sp_", sp_part)) {
      res <- gsub("_sp_", " sp. ", sp_part)
      res <- gsub("_", " ", res)
    } else if (grepl("SGB", sp_part)) {
      res <- gsub("_", " ", sp_part)
    } else {
      sub_parts <- strsplit(sp_part, "_")[[1]]
      if (length(sub_parts) >= 2) {
        genus <- sub_parts[1]
        species <- paste(sub_parts[-1], collapse = " ")
        genus_abbr <- paste0(substr(genus, 1, 1), ".")
        res <- paste(genus_abbr, species)
      } else {
        res <- gsub("_", " ", sp_part)
      }
    }

    if (sgb_id != "") paste(res, sgb_id) else res
  }, USE.NAMES = FALSE)
}

strength_label <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x >= 0.4 ~ "strong",
    x >= 0.2 ~ "moderate",
    TRUE ~ "weak"
  )
}

shape_group_from_overlap <- function(is_ob_da, overlap_class) {
  case_when(
    !is_ob_da ~ "Other_union_DA",
    overlap_class == "Ob-specific" ~ "Ob_specific",
    overlap_class == "Overlap with both Nw and Ow" ~ "BMI_overlap_all",
    overlap_class %in% c("Overlap with Nw", "Overlap with Ow") ~ "BMI_overlap_partial",
    TRUE ~ "Other_obese_DA"
  )
}

da_long <- bind_rows(lapply(names(feature_files), function(bmi_i) {
  read_csv(feature_files[[bmi_i]], show_col_types = FALSE) %>%
    mutate(source_bmi = bmi_i)
}))

ob_meta <- read_csv(feature_files[["Ob"]], show_col_types = FALSE) %>%
  transmute(
    feature,
    ob_clean_name = clean_name,
    Ob_DA_class = overlap,
    ob_coef_crc = coef_crc
  ) %>%
  distinct(feature, .keep_all = TRUE)

union_meta <- da_long %>%
  group_by(feature) %>%
  summarise(
    clean_name = dplyr::first(na.omit(clean_name)),
    source_bmi = paste(sort(unique(source_bmi)), collapse = ";"),
    coef_crc_pref = {
      vals <- coef_crc[!is.na(coef_crc)]
      if (length(vals) == 0) NA_real_ else vals[1]
    },
    .groups = "drop"
  ) %>%
  mutate(
    Union_DA_class = case_when(
      source_bmi == "Normal" ~ "Nw-spec",
      source_bmi == "Ow" ~ "Ow-spec",
      source_bmi == "Ob" ~ "Ob-spec",
      source_bmi == "Normal;Ow" ~ "Nw+Ow",
      source_bmi == "Normal;Ob" ~ "Nw+Ob",
      source_bmi == "Ob;Ow" ~ "Ow+Ob",
      source_bmi == "Normal;Ob;Ow" ~ "BMI-ind",
      TRUE ~ source_bmi
    ),
    crc_direction = if_else(coef_crc_pref >= 0, "CRC-enriched", "CRC-depleted")
  ) %>%
  left_join(ob_meta, by = "feature") %>%
  mutate(
    clean_name = coalesce(clean_name, ob_clean_name, format_sgb_names(feature)),
    is_obese_da = !is.na(Ob_DA_class),
    node_color_group = if_else(crc_direction == "CRC-enriched", "CRC_enriched", "CRC_depleted"),
    node_shape_group = shape_group_from_overlap(is_obese_da, Ob_DA_class)
  )

edges <- read_csv(edge_file, show_col_types = FALSE)

gain_classes <- c(
  "Positive correlation gained",
  "Negative correlation gained",
  "Positive correlation strengthened",
  "Negative correlation strengthened"
)

loss_classes <- c(
  "Positive correlation lost",
  "Negative correlation lost",
  "Positive correlation weakened",
  "Negative correlation weakened"
)

edge_sets <- list(
  Full = edges,
  Gain = edges %>% filter(edge_class %in% gain_classes),
  Loss = edges %>% filter(edge_class %in% loss_classes)
)

make_edges_long <- function(edge_df, set_name) {
  if (nrow(edge_df) == 0) return(tibble())

  source_annot <- union_meta %>%
    transmute(
      Gene1 = feature,
      source = clean_name,
      source_feature = feature,
      source_crc_direction = crc_direction,
      source_Ob_DA_class = if_else(is.na(Ob_DA_class), "Other union-DA", Ob_DA_class),
      source_Union_DA_class = Union_DA_class
    )

  target_annot <- union_meta %>%
    transmute(
      Gene2 = feature,
      target = clean_name,
      target_feature = feature,
      target_crc_direction = crc_direction,
      target_Ob_DA_class = if_else(is.na(Ob_DA_class), "Other union-DA", Ob_DA_class),
      target_Union_DA_class = Union_DA_class
    )

  base <- edge_df %>%
    left_join(source_annot, by = "Gene1") %>%
    left_join(target_annot, by = "Gene2") %>%
    mutate(
      candidate_set = paste0("Union_DAonly_ObeseDA_", set_name, "_ExcludingNwSameTrajectory"),
      same_nw_trajectory_removed = TRUE,
      Ob_specific_both_endpoints = da_side == "Both taxa are obese-DA",
      Ob_change_strength = case_when(
        edge_class %in% c("Positive correlation strengthened", "Negative correlation strengthened") ~ "Strengthened in CRC",
        edge_class %in% c("Positive correlation weakened", "Negative correlation weakened") ~ "Weakened in CRC",
        edge_class %in% c("Positive correlation gained", "Negative correlation gained") ~ "Gained in trajectory",
        edge_class %in% c("Positive correlation lost", "Negative correlation lost") ~ "Lost in trajectory",
        TRUE ~ edge_class
      )
    )

  stage_map <- tibble(
    Stage = c("CTRL", "AA", "CRC"),
    state_col = c("state__CTRL__Ob", "state__AA__Ob", "state__CRC__Ob"),
    cor_col = c("cor__CTRL__Ob", "cor__AA__Ob", "cor__CRC__Ob"),
    p_col = c("p_value__CTRL__Ob", "p_value__AA__Ob", "p_value__CRC__Ob")
  )

  bind_rows(lapply(seq_len(nrow(stage_map)), function(i) {
    st <- stage_map$Stage[i]
    state_col <- stage_map$state_col[i]
    cor_col <- stage_map$cor_col[i]
    p_col <- stage_map$p_col[i]

    base %>%
      filter(.data[[state_col]] != "0") %>%
      mutate(
        normal_ctrl_crc_diff_sig_p01 = !is.na(pValDiff_Normal) & pValDiff_Normal < 0.01,
        normal_ctrl_crc_diff_sig_p05 = !is.na(pValDiff_Normal) & pValDiff_Normal < 0.05,
        nw_progression_relation = case_when(
          !is.na(pValDiff_Normal) & pValDiff_Normal < 0.01 & trajectory__Normal != trajectory__Ob ~ "Shared progression, distinct trajectory",
          !is.na(pValDiff_Normal) & pValDiff_Normal < 0.01 & trajectory__Normal == trajectory__Ob ~ "Shared progression, same trajectory",
          TRUE ~ "Ob-specific progression"
        )
      ) %>%
      transmute(
        edge_id,
        source,
        target,
        source_feature,
        target_feature,
        source_crc_direction,
        target_crc_direction,
        source_Ob_DA_class,
        target_Ob_DA_class,
        source_Union_DA_class,
        target_Union_DA_class,
        candidate_set,
        Ob_specific_both_endpoints,
        Ob_change_strength,
        BMI = "Ob",
        trajectory = trajectory__Ob,
        trajectory_Normal = trajectory__Normal,
        state_CTRL = state__CTRL__Ob,
        state_AA = state__AA__Ob,
        state_CRC = state__CRC__Ob,
        Stage = st,
        rho = .data[[cor_col]],
        p_value = .data[[p_col]],
        edge_sign = if_else(.data[[cor_col]] >= 0, "positive", "negative"),
        abs_rho = abs(.data[[cor_col]]),
        edge_strength = strength_label(abs(.data[[cor_col]])),
        interaction = paste0(if_else(.data[[cor_col]] >= 0, "positive", "negative"), "_", st),
        Cytoscape_edge_name = paste(source, target, "Ob", st, sep = " | "),
        ctrl_crc_pValDiff = pValDiff,
        ctrl_crc_pValDiff_adj = pValDiff_adj,
        ctrl_crc_zScoreDiff = zScoreDiff,
        ctrl_crc_diff_sig_p001 = pValDiff < 0.001,
        ctrl_crc_diff_sig_p01 = pValDiff < 0.01,
        normal_ctrl_crc_pValDiff = pValDiff_Normal,
        normal_ctrl_crc_pValDiff_adj = pValDiff_adj_Normal,
        normal_ctrl_crc_diff_sig_p01,
        normal_ctrl_crc_diff_sig_p05,
        normal_edge_class = edge_class_normal,
        nw_progression_relation,
        edge_class,
        da_side
      )
  }))
}

make_nodes <- function(edge_df, set_name) {
  if (nrow(edge_df) == 0) return(tibble())

  nodes <- sort(unique(c(edge_df$Gene1, edge_df$Gene2)))

  degree_from_state <- function(node, state_cols) {
    edge_df %>%
      filter(Gene1 == node | Gene2 == node) %>%
      filter(if_any(all_of(state_cols), ~ .x != "0")) %>%
      summarise(n = n()) %>%
      pull(n)
  }

  degree_from_single <- function(node, state_col) {
    edge_df %>%
      filter(Gene1 == node | Gene2 == node) %>%
      filter(.data[[state_col]] != "0") %>%
      summarise(n = n()) %>%
      pull(n)
  }

  tibble(feature = nodes) %>%
    left_join(union_meta, by = "feature") %>%
    mutate(
      clean_name = coalesce(clean_name, format_sgb_names(feature)),
      node_id = clean_name,
      displayed_degree_total = map_int(feature, ~ sum(edge_df$Gene1 == .x | edge_df$Gene2 == .x)),
      displayed_degree_Nw = map_int(feature, ~ degree_from_state(.x, c("state__CTRL__Normal", "state__AA__Normal", "state__CRC__Normal"))),
      displayed_degree_Ow = map_int(feature, ~ degree_from_state(.x, c("state__CTRL__Ow", "state__AA__Ow", "state__CRC__Ow"))),
      displayed_degree_Ob = map_int(feature, ~ degree_from_state(.x, c("state__CTRL__Ob", "state__AA__Ob", "state__CRC__Ob"))),
      displayed_degree_CTRL = map_int(feature, ~ degree_from_single(.x, "state__CTRL__Ob")),
      displayed_degree_AA = map_int(feature, ~ degree_from_single(.x, "state__AA__Ob")),
      displayed_degree_CRC = map_int(feature, ~ degree_from_single(.x, "state__CRC__Ob")),
      candidate_set = paste0("Union_DAonly_ObeseDA_", set_name, "_ExcludingNwSameTrajectory")
    ) %>%
    transmute(
      node_id,
      feature,
      clean_name,
      crc_direction,
      Ob_DA_class = if_else(is.na(Ob_DA_class), "Other union-DA", Ob_DA_class),
      Union_DA_class,
      displayed_degree_total,
      displayed_degree_Nw,
      displayed_degree_Ow,
      displayed_degree_Ob,
      displayed_degree_CTRL,
      displayed_degree_AA,
      displayed_degree_CRC,
      node_color_group,
      node_shape_group,
      candidate_set
    ) %>%
    arrange(desc(displayed_degree_total), clean_name)
}

summary_rows <- list()
nwsum_rows <- list()

for (set_name in names(edge_sets)) {
  edge_df <- edge_sets[[set_name]]
  edges_long <- make_edges_long(edge_df, set_name)
  nodes_df <- make_nodes(edge_df, set_name)

  write_csv(edges_long, file.path(out_dir, paste0("Cytoscape_", set_name, "_Edges_Long.csv")))
  write_csv(nodes_df, file.path(out_dir, paste0("Cytoscape_", set_name, "_Nodes.csv")))

  summary_rows[[set_name]] <- tibble(
    Set = set_name,
    Edges = nrow(edge_df),
    Nodes = nrow(nodes_df),
    CTRL_Ob_edges = sum(edge_df$state__CTRL__Ob != "0"),
    AA_Ob_edges = sum(edge_df$state__AA__Ob != "0"),
    CRC_Ob_edges = sum(edge_df$state__CRC__Ob != "0")
  )

  nwsum_rows[[set_name]] <- edges_long %>%
    distinct(edge_id, nw_progression_relation, normal_edge_class) %>%
    count(nw_progression_relation, normal_edge_class, name = "n_edges") %>%
    mutate(Set = set_name, .before = 1)
}

write_csv(bind_rows(summary_rows), file.path(out_dir, "Cytoscape_Full_Gain_Loss_Summary.csv"))
write_csv(bind_rows(nwsum_rows), file.path(out_dir, "Cytoscape_Full_Gain_Loss_NwProgression_Summary.csv"))

cat("Saved Cytoscape-style exports to:\n", out_dir, "\n")
