## ============================================================
## Script 10: 5K → 10K WGCNA module gene redistribution figure
## ============================================================
## Generates a two-panel figure showing how each 5K WGCNA module
## (APOE4 microglia, 5,000 genes) maps to 10K modules (power=12).
##
## Panel A: Stacked horizontal bar chart (one bar per 5K module)
## Panel B: Gene-count heatmap table with row totals
##
## Inputs:
##   /mnt/shared-workspace/wgcna_gene_modules_noE2.csv
##   /mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv
##
## Output:
##   /mnt/results/WGCNA/WGCNA_10K/figures/module_split_5K_to_10K.png
## ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ── 1. Load data ──────────────────────────────────────────────────────────────
genes5k  <- read.csv("/mnt/shared-workspace/wgcna_gene_modules_noE2.csv",
                     stringsAsFactors = FALSE)
genes10k <- read.csv("/mnt/results/WGCNA/WGCNA_10K/tables/module_gene_assignments_10K_p12.csv",
                     stringsAsFactors = FALSE)

# ── 2. Cross-tabulation ───────────────────────────────────────────────────────
# Exclude 5K grey (not a real module)
genes5k_nonGrey <- genes5k %>% filter(module != "grey")

# Join to 10K assignments
joined <- genes5k_nonGrey %>%
  left_join(genes10k %>% select(ensembl_id, module_10k = module), by = "ensembl_id") %>%
  mutate(dest = case_when(
    is.na(module_10k)       ~ "Lost",   # not in 10K dataset
    module_10k == "grey"    ~ "Lost",   # assigned to 10K grey
    TRUE                    ~ module_10k
  ))

long_data <- joined %>%
  count(module_5k = module, dest) %>%
  group_by(module_5k) %>%
  mutate(total = sum(n)) %>%
  ungroup()

# ── 3. Ordering ───────────────────────────────────────────────────────────────
module_order <- long_data %>%
  distinct(module_5k, total) %>%
  arrange(desc(total)) %>%
  pull(module_5k)

dest_order_base <- long_data %>%
  filter(dest != "Lost") %>%
  group_by(dest) %>%
  summarise(total_dest = sum(n)) %>%
  arrange(desc(total_dest)) %>%
  pull(dest)

dest_order <- c(dest_order_base, "Lost")

# Destinations to show in Panel B (>= 5 genes from any 5K module)
keep_dests <- long_data %>%
  group_by(dest) %>%
  summarise(max_n = max(n)) %>%
  filter(max_n >= 5 | dest == "Lost") %>%
  pull(dest)

# ── 4. Colour mapping ─────────────────────────────────────────────────────────
module_colours <- c(
  "turquoise"   = "#00CED1",
  "blue"        = "#4169E1",
  "brown"       = "#8B4513",
  "yellow"      = "#FFD700",
  "green"       = "#228B22",
  "red"         = "#DC143C",
  "black"       = "#1A1A1A",
  "pink"        = "#FF69B4",
  "magenta"     = "#CC00CC",
  "purple"      = "#7B2D8B",
  "greenyellow" = "#ADFF2F",
  "salmon"      = "#FA8072",
  "tan"         = "#D2B48C",
  "cyan"        = "#00BFFF",
  "Lost"        = "#AAAAAA"
)

# ── 5. Panel A: stacked horizontal bar chart ──────────────────────────────────
long_data_A <- long_data %>%
  mutate(
    module_5k = factor(module_5k, levels = rev(module_order)),
    dest      = factor(dest, levels = dest_order),
    label     = ifelse(n >= 15, as.character(n), "")
  )

panelA <- ggplot(long_data_A, aes(x = n, y = module_5k, fill = dest)) +
  geom_col(position = "stack", width = 0.7, colour = "white", linewidth = 0.3) +
  geom_col(
    data = long_data_A %>% filter(module_5k == "yellow"),
    aes(x = n, y = module_5k),
    position = "stack", width = 0.7,
    fill = NA, colour = "black", linewidth = 0.8, inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.8, colour = "white", fontface = "bold"
  ) +
  scale_fill_manual(values = module_colours, name = "10K module destination",
                    guide = guide_legend(ncol = 1)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1400, 200)) +
  labs(
    x = "Number of genes", y = "5K module",
    title = "A   5K → 10K WGCNA module gene redistribution",
    subtitle = paste0(
      "Bar length = total genes in 5K module. Fill = 10K destination. ",
      "'Lost' = 10K grey or absent from 10K dataset.\n",
      "Yellow bar (biological focus) outlined in black."
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(
      colour = c("magenta3","hotpink","red3","darkgreen","brown",
                 "gold3","royalblue","saddlebrown","turquoise4"),
      face = "bold", size = 11
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
  )

# ── 6. Panel B: gene-count heatmap table ──────────────────────────────────────
real_dests <- setdiff(as.character(keep_dests), "Lost")
all_cols   <- c(real_dests, "Lost", "Total")

row_totals_df <- long_data %>%
  group_by(module_5k) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(dest = "Total", fill_val = NA_real_, cell_label = as.character(n),
         is_yellow = module_5k == "yellow")

table_data <- long_data %>%
  filter(dest %in% keep_dests) %>%
  group_by(module_5k) %>%
  mutate(fill_val = ifelse(dest == "Lost", NA_real_,
                           n / max(n[dest != "Lost"], na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(cell_label = ifelse(n == 0, "", as.character(n)),
         is_yellow  = module_5k == "yellow")

table_data <- bind_rows(table_data, row_totals_df) %>%
  mutate(
    dest      = factor(dest, levels = all_cols),
    module_5k = factor(module_5k, levels = rev(module_order))
  )

x_colours <- c(
  "turquoise" = "turquoise4", "brown" = "saddlebrown", "blue" = "royalblue",
  "green" = "darkgreen", "pink" = "hotpink", "purple" = "purple4",
  "yellow" = "gold3", "greenyellow" = "chartreuse4", "salmon" = "salmon3",
  "magenta" = "magenta3", "black" = "black", "red" = "red3",
  "tan" = "tan4", "cyan" = "deepskyblue3", "Lost" = "grey40", "Total" = "grey20"
)
x_col_vec <- x_colours[all_cols]

panelB <- ggplot(table_data, aes(x = dest, y = module_5k)) +
  geom_tile(
    data = table_data %>% filter(!dest %in% c("Lost", "Total")),
    aes(fill = fill_val), colour = "white", linewidth = 0.5
  ) +
  geom_tile(
    data = table_data %>% filter(dest == "Lost"),
    fill = "#DDDDDD", colour = "white", linewidth = 0.5
  ) +
  geom_tile(
    data = table_data %>% filter(dest == "Total"),
    fill = "#FFF8DC", colour = "grey70", linewidth = 0.5
  ) +
  geom_tile(
    data = table_data %>% filter(module_5k == "yellow"),
    colour = "black", linewidth = 1.0, fill = NA
  ) +
  geom_vline(xintercept = length(real_dests) + 0.5,
             colour = "grey50", linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = length(real_dests) + 1.5,
             colour = "grey30", linewidth = 0.8) +
  geom_text(aes(label = cell_label), size = 3.0, colour = "grey10", fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#2166AC", na.value = "#DDDDDD",
                      name = "Relative count\n(row-normalised,\nexcl. Lost)") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "B   Gene counts per 5K → 10K module destination") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0, size = 9, face = "bold",
                                colour = x_col_vec),
    axis.text.y = element_text(
      colour = c("magenta3","hotpink","red3","darkgreen","brown",
                 "gold3","royalblue","saddlebrown","turquoise4"),
      face = "bold", size = 10
    ),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 11),
    panel.grid = element_blank()
  )

# ── 7. Combine and save ───────────────────────────────────────────────────────
combined <- panelA / panelB +
  plot_layout(heights = c(3, 2)) +
  plot_annotation(
    title    = "5K → 10K WGCNA module gene redistribution",
    subtitle = "APOE4 microglia dataset | 5K run: 9 modules (grey excluded) | 10K run: power=12, 14 modules",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

out_path <- "/mnt/results/WGCNA/WGCNA_10K/figures/module_split_5K_to_10K.png"
ggsave(out_path, combined, width = 15, height = 11, dpi = 150, bg = "white")
message("Saved: ", out_path)
