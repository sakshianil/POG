#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
library(tidyverse)
library(patchwork)

option_list <- list(
  make_option("--input", type="character", default=NULL, help="Input CSV file"),
  make_option("--output", type="character", default="results/aminoacid_composition", help="Output directory")
)

opt <- parse_args(OptionParser(option_list=option_list))
input_file <- opt$input
output_dir <- opt$output
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read input
df <- read.csv(input_file, sep = ";", stringsAsFactors = FALSE)

# Process regions function
process_region <- function(data, region_col, label) {
    data %>%
        filter(!is.na(.data[[region_col]])) %>%
        select(GeneID, SecretionStatus, all_of(region_col)) %>%
        mutate(sequence = .data[[region_col]]) %>%
        select(-all_of(region_col)) %>%
        mutate(sequence = strsplit(sequence, "")) %>%
        unnest_longer(sequence, values_to = "AA") %>%
        group_by(GeneID, SecretionStatus) %>%
        mutate(Position = row_number(), Region = label) %>%
        ungroup() %>%
        filter(AA %in% c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
}
# Process signal peptide and conserved region
sp <- process_region(df, "SignalPeptide", "Signal Peptide")
cr <- process_region(df, "ConservedRegion", "Conserved Region") %>%
    mutate(Position = Position + 30)

all_data <- bind_rows(sp, cr)

aa_groups <- tribble(
    ~AA, ~Group,
    "K", "Positively charged\n[K/R/H]",
    "R", "Positively charged\n[K/R/H]",
    "H", "Positively charged\n[K/R/H]",
    "D", "Negatively charged\n[D/E]",
    "E", "Negatively charged\n[D/E]",
    "S", "Polar uncharged\n[S/T/N/Q]",
    "T", "Polar uncharged\n[S/T/N/Q]",
    "N", "Polar uncharged\n[S/T/N/Q]",
    "Q", "Polar uncharged\n[S/T/N/Q]",
    "A", "Small/other\n[A/G/P/C/M]",
    "G", "Small/other\n[A/G/P/C/M]",
    "P", "Small/other\n[A/G/P/C/M]",
    "C", "Small/other\n[A/G/P/C/M]",
    "M", "Small/other\n[A/G/P/C/M]",
    "Y", "Aromatic\n[Y/F/W]",
    "F", "Aromatic\n[Y/F/W]",
    "W", "Aromatic\n[Y/F/W]",
    "L", "Hydrophobic\n[L/I/V]",
    "I", "Hydrophobic\n[L/I/V]",
    "V", "Hydrophobic\n[L/I/V]"
)

# Prepare heatmap data
heatmap_df <- all_data %>%
    group_by(SecretionStatus, Region, Position, AA) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(SecretionStatus, Region, Position) %>%
    mutate(Frequency = Count / sum(Count)) %>%
    ungroup()

heatmap_df <- heatmap_df %>% left_join(aa_groups, by = "AA")

heatmap_df$Region <- factor(heatmap_df$Region, levels = c("Signal Peptide", "Conserved Region"))
heatmap_df$SecretionStatus <- factor(heatmap_df$SecretionStatus, levels = c("Secreted", "Non-Secreted", "Validated Effectors"))
heatmap_df$AA <- factor(heatmap_df$AA, levels = c("K", "R", "H", "D", "E", "S", "T", "N", "Q", "A", "G", "P", "C", "M", "Y", "F", "W", "L", "I", "V"))
# Calculate percentages for bar plot before merging top AA
group_summary <- heatmap_df %>%
    group_by(SecretionStatus, Group, Region) %>%
    summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(SecretionStatus, Region) %>%
    mutate(Percentage = Count / sum(Count) * 100) %>%
    ungroup()

# Transform frequency for plotting
heatmap_df <- heatmap_df %>%
    mutate(Frequency_transformed = sqrt(Frequency))
# Calculate percentages for amino acids across all positions per SecretionStatus
aa_percentage_labels <- heatmap_df %>%
    group_by(SecretionStatus, AA) %>%
    summarise(Total_Frequency = sum(Frequency), .groups = "drop") %>%
    group_by(SecretionStatus) %>%
    mutate(Percentage = Total_Frequency / sum(Total_Frequency) * 100) %>%
    mutate(Label = paste0(round(Percentage, 1), ""))
# Merge labels back into heatmap_df
heatmap_df <- heatmap_df %>%
    left_join(aa_percentage_labels %>% select(SecretionStatus, AA, Label), by = c("SecretionStatus", "AA")) %>%
    mutate(Frequency_transformed = Frequency ^ 0.3)


p <- ggplot(heatmap_df, aes(x = Position, y = AA, fill = Frequency_transformed)) +
    geom_tile(color = NA, na.rm = TRUE) +

    scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#D7191C",
                         midpoint = 0.15, name = "AA Frequency\n(sqrt scaled)") +

    facet_wrap2(~ SecretionStatus, ncol = 1, strip.position = "right", scales = "fixed") +
    geom_text(data = aa_percentage_labels, aes(x = 155, y = AA, label = Label),
              inherit.aes = FALSE, size = 3) +
    theme_minimal(base_size = 10) +
  scale_y_discrete(expand = c(0, 0)) +
    theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text = element_text(size = 12, face = "bold", family = "Arial"),
        axis.text.x = element_text(size = 10, face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 8, angle= 10, face = "bold", family = "Arial"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        strip.placement = "outside",
        panel.spacing = unit(0.2, "lines"),
        legend.position = "right",
        plot.margin = margin(3, 3, 3, 3)
    ) +
    labs(x = "Amino Acid Position in Protein Sequence", y = "Amino Acid")



top_aa_per_group_region <- heatmap_df %>%
    group_by(SecretionStatus, Group, Region, AA) %>%
    summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    group_by(SecretionStatus, Group, Region) %>%
    slice_max(Total_Count, n = 1, with_ties = FALSE) %>%
    ungroup()
group_summary <- group_summary %>%
    left_join(top_aa_per_group_region %>% select(SecretionStatus, Group, Region, TopAA = AA),
              by = c("SecretionStatus", "Group", "Region"))

o_group_summary <- ggplot(group_summary, aes(x = Group, y = Percentage, fill = SecretionStatus)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = paste0(round(Percentage, 1), "\n(", TopAA, ")")),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 3, family = "Arial") +
    scale_fill_brewer(palette = "Set2") +
    facet_wrap(~Region, ncol = 1, scales = "free_y") +  # Separate SP and CR panels
    theme_minimal(base_size = 11, base_family = "Arial") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
        title = "Amino Acid Composition",
        y = "Percentage (%)",
        x = "Amino Acid Group"
    )

# Only draw SP + Conserved bar once using patchwork

sp_cons_bar <- ggplot() +
    annotate("rect", xmin = 0, xmax = 30, ymin = 0, ymax = 1, fill = "white", color = "black") +
    annotate("rect", xmin = 31, xmax = 150, ymin = 0, ymax = 1, fill = "white", color = "black") +
    annotate("text", x = 15, y = 0.5, label = "SP", size = 4, fontface = "bold", family = "Arial") +
    annotate("text", x = 90, y = 0.5, label = "Conserved region", size = 4, fontface = "bold", family = "Arial") +
    theme_void() +
    coord_cartesian(xlim = c(0, 150), ylim = c(0, 1), expand = FALSE)

# Combine using patchwork: SP + Conserved on top, heatmap below
g <- sp_cons_bar / p + plot_layout(heights = c(0.05, 1)) +
    annotate("segment", x = 30, xend = 30, y = -Inf, yend = Inf,
             color = "black", linetype = "dashed", linewidth = 0.3) +
    annotate("segment", x = 31, xend = 31, y = -Inf, yend = Inf,
             color = "black", linetype = "dashed", linewidth = 0.3)

final_layout_extended <- ((a)/(e | b) /(c | d) | f | (g)) +
    plot_layout(guides = "collect", heights = c(1.5, 1.2, 1, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        plot.margin = margin(2, 2, 2, 2)
)

ggsave(file.path(output_dir, "fig_aa_comp.tiff"), 
       final_layout_extended, width = 20, height = 25, 
       units = "cm", dpi = 600, compression = "lzw")



############Supp figure for all species seperately generated for aminoa acid composition calculation#######

species_folders <- c("harab_up", "pult_up", "phal_up", "pinf_up", "psoj_up")
species_names <- c("Harab", "Pult", "Phal", "Pinf", "Psoj")

# Read CSVs and add species label
csv_list <- list()
for (i in seq_along(species_folders)) {
    path <- file.path(base_path, species_folders[i], paste0(strsplit(species_folders[i], "_")[[1]][1], "_sp_cr.csv"))
    if (file.exists(path)) {
        df <- read.csv(path, stringsAsFactors = FALSE)
        df$SpeciesStatus <- paste0(species_names[i], "_Secreted")
        csv_list[[i]] <- df
    }
}

# Add validated effectors
validated_path <- file.path(base_path, "validated_effectors_secreted_non_secrted", "protein_regions_all_types.csv")
validated_df <- read.csv(validated_path, sep = ";", stringsAsFactors = FALSE)
validated_df$SpeciesStatus <- "Validated_Effectors"

# Combine everything
all_df <- bind_rows(csv_list, validated_df)
# Process function
process_region <- function(data, region_col, label) {
    data %>%
        filter(!is.na(.data[[region_col]])) %>%
        mutate(sequence = strsplit(.data[[region_col]], "")) %>%
        unnest_longer(sequence, values_to = "AA") %>%
        group_by(ProteinID = coalesce(GeneID, ProteinID), SpeciesStatus) %>%
        mutate(Position = row_number(), Region = label) %>%
        ungroup() %>%
        filter(AA %in% c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                         "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"))
}

sp <- process_region(all_df, "SignalPeptide", "Signal Peptide")
cr <- process_region(all_df, "ConservedRegion", "Conserved Region") %>%
    mutate(Position = Position + 30)

all_data <- bind_rows(sp, cr)
# Set factor levels for species in the desired order and with pretty names
species_levels <- c("Harab_Secreted", "Phal_Secreted", "Pinf_Secreted", "Psoj_Secreted", "Pult_Secreted", "Validated_Effectors")
species_labels <- c("*H. arabidopsidis*", "*Pl. halstedii*", "*P. infestans*", "*P. sojae*", "*P. ultimum*", "*Validated Effectors*")
# Prepare heatmap data
heatmap_df <- all_data %>%
    group_by(SpeciesStatus, Region, Position, AA) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(SpeciesStatus, Region, Position) %>%
    mutate(Frequency = Count / sum(Count)) %>%
    ungroup()

# Calculate AA percentages for labels
aa_percentage_labels <- heatmap_df %>%
    group_by(SpeciesStatus, AA) %>%
    summarise(Total_Frequency = sum(Frequency), .groups = "drop") %>%
    group_by(SpeciesStatus) %>%
    mutate(Percentage = Total_Frequency / sum(Total_Frequency) * 100) %>%
    mutate(Label = paste0(round(Percentage, 1), "%"))

# Plotting using facet_wrap
library(ggplot2)


p <- ggplot(heatmap_df, aes(x = Position, y = AA, fill = Frequency ^ 0.3)) +
    geom_tile(color = NA, na.rm = TRUE) +
    scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#D7191C",
                         midpoint = 0.15, name = "AA Frequency\n(sqrt scaled)") +
    facet_wrap(~SpeciesStatus, ncol = 1, strip.position = "right", scales = "fixed") +
    geom_text(data = aa_percentage_labels, aes(x = 155, y = AA, label = Label),
              inherit.aes = FALSE, size = 3) +
    theme_minimal(base_size = 10) +
    labs(x = "Amino Acid Position in Protein Sequence", y = "Amino Acid") +
    theme(strip.text.y.left = element_markdown(size = 12, face = "italic"),
          strip.placement = "outside",
      panel.spacing = unit(0.2, "lines"))

# Add SP + Conserved bar using patchwork layout
sp_cons_bar <- ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = 30, ymin = 0, ymax = 1, fill = "white", color = "black") +
    ggplot2::annotate("rect", xmin = 31, xmax = 150, ymin = 0, ymax = 1, fill = "white", color = "black") +
    ggplot2::annotate("text", x = 15, y = 0.5, label = "SP", size = 4, fontface = "bold") +
    ggplot2::annotate("text", x = 90, y = 0.5, label = "Conserved region", size = 4, fontface = "bold") +
    theme_void()

combined_plot <- sp_cons_bar / p + plot_layout(heights = c(0.05, 1))

ggsave(file.path(output_dir, "fig_supp_aacomp.tiff"), combined_plot, width = 20, height = 25, units = "cm", dpi = 600, compression = "lzw")
