library(tidyverse)
library(ggrepel)

setwd("~/Desktop/Stomentosus_PopGen/PCA/")

# load data from plink
pca <- read_table("pca_results.eigenvec", col_names = FALSE)
eigenval <- scan("pca_results.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


# remake data.frame
pca <- as_tibble(data.frame(pca))

# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)


# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_minimal()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(PC1, PC2)) + geom_point(size = 3)
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


# add the isolate information
info <- read.csv("~/Desktop/Stomentosus_PopGen/popgen_sample_metadata.csv", header = T)
info <- na.omit(info)
pca2 <- merge(pca, info, by="ind", all = TRUE)


c <- ggplot(pca2, aes(PC1, PC2, col = site)) + geom_point(size = 5)

c <- c + coord_equal() + 
  theme_bw() + 
  theme(axis.text  = element_text(size = 16),
                              axis.title.y = element_text(size = 16),
                              axis.title.x = element_text(size = 16),
                              plot.background = element_blank(), 
                              panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), 
                              legend.position = "none") 

c <- c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

c 

# highlight the isolates in the zn20 group (top 10 and bottom 10 Zn tol)
# and mark them with open circles
# also have an arrow pointing at ZnT (230) and ZnS (220)

pca2

pca2 %>%
  slice_min(zn_ec_mm, n = 10)

pca2 %>%
  slice_max(zn_ec_mm, n = 10)



pc1_col <- "PC1"         # change if your PC1 column is named differently
pc2_col <- "PC2"         # change if your PC2 column is named differently
id_col  <- "ind"     # the column with isolate names/IDs
metric_col <- "zn_ec_mm"   # the column with the numeric metric to rank by

top_n <- 10              # how many to highlight at the top
bottom_n <- 10           # how many to highlight at the bottom
labels_to_annotate <- tibble(
  ind = c("ACCCACAACT-CGGTGGTAAG", "GTATGGCTTC-CCGATCGCCA"),
  label = c("ZnS", "ZnT")
)  # <-- replace with the two names you want to label

# ---------------------------
# BUILD TOP/BOTTOM SETS
# ---------------------------
top_ids <- pca2 %>%
  slice_max(order_by = !!sym(metric_col), n = top_n, with_ties = FALSE) %>%
  pull(!!sym(id_col))

bottom_ids <- pca2 %>%
  slice_min(order_by = !!sym(metric_col), n = bottom_n, with_ties = FALSE) %>%
  pull(!!sym(id_col))

# Add a category for highlighting
pca2 <- pca2 %>%
  mutate(
    highlight = case_when(
      .data[[id_col]] %in% top_ids    ~ "Top 10",
      .data[[id_col]] %in% bottom_ids ~ "Bottom 10",
      TRUE                            ~ "None"
    )
  )

# Subset for labels
label_df <- pca2 %>% filter(.data[[id_col]] %in% labels_to_annotate$ind)

library(grid)  # for arrow()

# --- Identify your two points and set label names ---
label_df_ZnS <- pca2 %>%
  filter(ind == "ACCCACAACT-CGGTGGTAAG") %>%
  mutate(label_name = "ZnS")

label_df_ZnT <- pca2 %>%
  filter(ind == "GTATGGCTTC-CCGATCGCCA") %>%
  mutate(label_name = "ZnT")

# --- Plot ranges for scale-aware nudges ---
x_rng <- diff(range(pca2$PC1, na.rm = TRUE))
y_rng <- diff(range(pca2$PC2, na.rm = TRUE))

# Position labels farther away to get long arrows
nudge_x_S <- 0.18  # left
nudge_y_S <- 0.14  # down
nudge_x_T <- 0.22  # right
nudge_y_T <- 0.12  # down

label_df_ZnS <- label_df_ZnS %>%
  mutate(x_lab = PC1 - nudge_x_S * x_rng,
         y_lab = PC2 - nudge_y_S * y_rng)

label_df_ZnT <- label_df_ZnT %>%
  mutate(x_lab = PC1 + nudge_x_T * x_rng,
         y_lab = PC2 - nudge_y_T * y_rng)

# --- Add a small gap between the label and the arrow segment start ---
# gap is a fraction of the overall plot range (adjust to taste)
gap_frac <- 0.02  # ~2% of the diagonal range

# Helper: compute a small offset from the label towards the point
add_gap <- function(df, gap_frac, x_rng, y_rng) {
  # vector from label to point
  dx <- df$PC1 - df$x_lab
  dy <- df$PC2 - df$y_lab
  # normalize to unit length in plot coordinates
  len <- sqrt((dx / x_rng)^2 + (dy / y_rng)^2)
  # avoid division by zero
  len[len == 0] <- 1
  # move start point a bit towards the point
  df$x_start <- df$x_lab + (dx / len) * gap_frac * x_rng
  df$y_start <- df$y_lab + (dy / len) * gap_frac * y_rng
  df
}

label_df_ZnS <- add_gap(label_df_ZnS, gap_frac, x_rng, y_rng)
label_df_ZnT <- add_gap(label_df_ZnT, gap_frac, x_rng, y_rng)

# --- Build plot ---
c <- ggplot(pca2, aes(PC1, PC2, color = site)) +
  # Base points colored by site
  geom_point(size = 5, alpha = 0.9) +
  # Open-circle overlays for Top/Bottom 10
  geom_point(
    data = subset(pca2, highlight != "None"),
    shape = 21, fill = NA, size = 6.2, stroke = 1.2
  ) +
  # Text labels at nudged positions
  geom_text(
    data = label_df_ZnS,
    aes(x = x_lab, y = y_lab, label = label_name),
    color = "#E69F00", size = 5, fontface = "bold"
  ) +
  geom_text(
    data = label_df_ZnT,
    aes(x = x_lab, y = y_lab, label = label_name),
    color = "#0072B2", size = 5, fontface = "bold"
  ) +
  # Arrow segments starting slightly away from the text (to avoid overlap)
  geom_segment(
    data = label_df_ZnS,
    aes(x = x_start, y = y_start, xend = PC1, yend = PC2),
    color = "#E69F00", size = 0.6, lineend = "butt",
    arrow = arrow(length = unit(0.010, "npc"), type = "closed")
  ) +
  geom_segment(
    data = label_df_ZnT,
    aes(x = x_start, y = y_start, xend = PC1, yend = PC2),
    color = "#0072B2", size = 0.6, lineend = "butt",
    arrow = arrow(length = unit(0.010, "npc"), type = "closed")
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    axis.text  = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

c
