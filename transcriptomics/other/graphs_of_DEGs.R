library(ggplot2)
library(dplyr)
library(scales)

# read in rds from DEG analysis - or use object from 01_DEG script
# res_shrunken_df <- readRDS("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomtomento230_ref/Tissue/220deseq_results_shrunken_lfc_df.rds")

# here use res_df from the DEseq
deg_data <- res_shrunken_df %>%
  mutate(
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1.5 ~ "up",
      padj < 0.05 & log2FoldChange < -1.5 ~ "down",
      TRUE ~ "not_significant"
    )
  )

# for tissue exclude 5 mM - but not for bioassay
#deg_data <- deg_data[!deg_data$result_name=="treatment_5_vs_0",]

# Set your desired x-axis order
deg_data$result_name <- factor(deg_data$result_name, 
                               levels = c(
                                 #"treatment_0.1_vs_0",
                                 "treatment_1_vs_0", 
                                 #"treatment_2.5_vs_0",
                                 "treatment_5_vs_0",
                                 "treatment_10_vs_0"))

# Count up/down regulated genes per condition
reg_counts <- deg_data %>%
  filter(regulation != "not_significant") %>%
  group_by(result_name, regulation) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(regulation == "down", -count, count))

# Scale counts to match y-axis range
ymin <- -10
ymax <- 10
max_count <- max(abs(reg_counts$count))
scale_factor <- ymax / max_count
reg_counts <- reg_counts %>%
  mutate(scaled_count = count * scale_factor)

# Separate significant and non-significant data
sig_data <- deg_data %>% filter(regulation != "not_significant")
nonsig_data <- deg_data %>% filter(regulation == "not_significant")

# Plot
p <- ggplot() +
  # Non-significant points
  geom_jitter(data = nonsig_data, aes(x = result_name, y = log2FoldChange, size = -log10(padj)),
              color = "grey", width = 0.2, alpha = 0.2) +
  # Significant points
  geom_jitter(data = sig_data, aes(x = result_name, y = log2FoldChange, color = regulation, size = -log10(padj)),
              width = 0.2, alpha = 0.8) +
  # DEG count bars
  geom_bar(data = reg_counts, aes(x = result_name, y = scaled_count, fill = regulation),
           stat = "identity", position = "identity", alpha = 0.2, width = 0.6) +
  # DEG count labels
  geom_text(
    data = reg_counts,
    aes(x = result_name, y = scaled_count, label = abs(count)),
    vjust = ifelse(reg_counts$count > 0, -1, 1.5),  # Above for up, below for down
    hjust = 0.5,  # Centered horizontally
    size = 4.5,
    inherit.aes = FALSE
  ) +
  # Scales and themes
  # "220" = "#E69F00", "230" = "#0072B2"
  scale_color_manual(values = c("up" = "#E69F00", "down" = "#E69F00")) +
  scale_fill_manual(values = c("up" = "#D55E00", "down" = "#56B4E9")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    limits = c(ymin, ymax),
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of Significant Genes")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(
    x = "Zn Treatment (mM)",
    y = expression(Log[2]~Fold~Change),
    color = "Regulation",
    size = "-log10(padj)",
    fill = "Regulation"
  )

p

# make sure all x values are present even if no data
# and ensure same p val scale across isoaltes
p <- p + scale_x_discrete(limits = c(
  #"treatment_0.1_vs_0", 
  "treatment_1_vs_0",
  #"treatment_2.5_vs_0",
  "treatment_5_vs_0",
  "treatment_10_vs_0")) +
  
  scale_size_continuous(
    name = "-log10(padj)",
    limits = c(0, 20),         # Set the same range across datasets
    breaks = c(5, 10, 15, 20), # Consistent legend ticks
    range = c(1, 5)            # Visual size range of points
  )


p + theme(
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 16),
    #axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.x = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title = element_blank(),
    #axis.title.x.bottom = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_blank()) +
  #scale_x_discrete(labels = c("treatment_0.1_vs_0" = "0.1",
                              "treatment_1_vs_0" = "1",
                              "treatment_5_vs_0" = "5",
                              "treatment_10_vs_0" = "10"))
                              


####
## updated version


library(ggplot2)
library(dplyr)
library(scales)

# --- 0) Precompute regulation & a stable size variable -----------------------
deg_data <- res_shrunken_df %>%
  mutate(
    regulation = case_when(
      padj < 0.05 & log2FoldChange >  1.5 ~ "up",
      padj < 0.05 & log2FoldChange < -1.5 ~ "down",
      TRUE ~ "not_significant"
    ),
    # Clean padj to avoid Inf and NA in -log10
    padj_clean = case_when(
      is.na(padj) ~ NA_real_,
      padj == 0  ~ 1e-300,           # floor to avoid Inf
      TRUE       ~ padj
    ),
    neglog10_p = -log10(padj_clean),
    # Cap at global max to ensure comparability across runs
    neglog10_p = pmin(neglog10_p, 20)   # <- choose your global cap
  )

# Optional: exclude certain contrasts (keep as you had it)
# deg_data <- deg_data[!deg_data$result_name == "treatment_5_vs_0", ]

# Fixed x-levels
deg_data$result_name <- factor(
  deg_data$result_name,
  levels = c("treatment_1_vs_0", "treatment_5_vs_0", "treatment_10_vs_0")
)

# --- 1) Count up/down DEGs per contrast and scale to y2 axis -----------------
reg_counts <- deg_data %>%
  filter(regulation != "not_significant") %>%
  group_by(result_name, regulation) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(count = ifelse(regulation == "down", -count, count))

ymin <- -10
ymax <-  10
max_count  <- max(abs(reg_counts$count), na.rm = TRUE)
scale_factor <- if (max_count > 0) ymax / max_count else 1

reg_counts <- reg_counts %>%
  mutate(scaled_count = count * scale_factor)

# --- 2) Split sig/nonsig for plotting layers --------------------------------
sig_data    <- filter(deg_data, regulation != "not_significant")
nonsig_data <- filter(deg_data, regulation == "not_significant")

# --- 3) Plot -----------------------------------------------------------------
p <- ggplot() +
  # Non-significant points (uniform color/alpha)
  geom_jitter(
    data = nonsig_data,
    aes(x = result_name, y = log2FoldChange, size = neglog10_p),
    color = "grey50", width = 0.2, alpha = 0.2, na.rm = TRUE
  ) +
  # Significant points (colored by regulation)
  geom_jitter(
    data = sig_data,
    aes(x = result_name, y = log2FoldChange, color = regulation, size = neglog10_p),
    width = 0.2, alpha = 0.8, na.rm = TRUE
  ) +
  # DEG count bars on secondary axis scale
  geom_bar(
    data = reg_counts,
    aes(x = result_name, y = scaled_count, fill = regulation),
    stat = "identity", position = "identity", alpha = 0.2, width = 0.6
  ) +
  # DEG count labels
  geom_text(
    data = reg_counts,
    aes(x = result_name, y = scaled_count, label = abs(count)),
    vjust = ifelse(reg_counts$count > 0, -1, 1.5),
    hjust = 0.5, size = 4.5, inherit.aes = FALSE
  ) +
  # Colors (you can tweak to your palette)
  # "220" = "#E69F00", "230" = "#0072B2"
  scale_color_manual(values = c("up" = "#0072B2", "down" = "#0072B2")) +
  scale_fill_manual(values  = c(up = "#D55E00", down = "#56B4E9")) +
  
  # --- The critical bit: a single, fixed size scale shared by all layers -----
scale_size_continuous(
  name   = "-log10(padj)",
  limits = c(0, 20),                  # fixed across runs
  breaks = c(5, 10, 15, 20),
  labels = c("5", "10", "15", "20"),
  range  = c(1, 5)) +           # visual size range (same across runs)
  # Y scales (primary + secondary)
  
  scale_y_continuous(
    limits = c(ymin, ymax),
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name   = "Number of Significant Genes",
      labels = function(b) scales::number(abs(b), accuracy = 1)
    )
  ) +
  
  
  # X scale: ensure all contrasts appear even if empty
  scale_x_discrete(limits = c("treatment_1_vs_0", "treatment_5_vs_0", "treatment_10_vs_0")) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text   = element_text(size = 10),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"  # set to "none" if you truly don’t want the size legend
  ) +
  labs(
    x = "Zn Treatment (mM)",
    y = expression(Log[2]~Fold~Change),
    color = "Regulation",
    fill  = "Regulation"
  )

p


