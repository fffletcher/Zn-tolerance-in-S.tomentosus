library(tidyverse)
library(qqman)

setwd("~/Desktop/Stomentosus_PopGen/FST/")

# read in the file
fst <- read_tsv("top_bottom.fst.weir.fst")
# Remove rows with NA
fst <- fst[complete.cases(fst), ]
# change all negative Fst values to 0
fst$WEIR_AND_COCKERHAM_FST[fst$WEIR_AND_COCKERHAM_FST < 0] <- 0
# Convert CHROM to numeric scaffold IDs
# Remove 'scaffold_' prefix and convert to numeric
fst$CHR <- as.numeric(gsub("scaffold_", "", fst$CHROM))
# Add SNP ID - useful for highlighting
fst$SNP <- 1:nrow(fst)


# windowed 5kb Fst
fst_windowed <- read_tsv("top_bottom.fst.windowed.weir.fst")

# Remove 'scaffold_' prefix and convert to numeric
fst_windowed$CHR <- as.numeric(gsub("scaffold_", "", fst_windowed$CHROM))

# Add SNP ID (can be row number)
fst_windowed$SNP <- paste0("win_", 1:nrow(fst_windowed))

# Use BIN_START as the position
fst_windowed$BP <- fst_windowed$BIN_START

# plot
manhattan(fst_windowed,
         chr = "CHR",
         bp = "BP",
         p = "WEIGHTED_FST",
         snp = "SNP",
         logp = FALSE,
         ylab = "Weighted FST",
         main = "Manhattan Plot of Windowed FST (5kb)")

# which windows have an fst of 1

fst_one <- fst_windowed %>%
  filter(WEIGHTED_FST == 1)
# how many windows have FST = 1
nrow(fst_one)


# Calculate the top 1% threshold including FST = 1
fst_threshold_1p <- quantile(fst_windowed$WEIGHTED_FST, 0.99)
fst_threshold_5p <- quantile(fst_windowed$WEIGHTED_FST, 0.95)

# Identify top % windows including FST = 1
fst_top_1 <- fst_windowed %>% filter(WEIGHTED_FST >= fst_threshold_1p)
fst_top_5 <- fst_windowed %>% filter(WEIGHTED_FST >= fst_threshold_5p)

# how many windows in the top 1%
nrow(fst_top_1)
# how many windows in the top 5%
nrow(fst_top_5)
# what is the lowest FST in these groups
min(fst_top_1$WEIGHTED_FST, na.rm = TRUE)
min(fst_top_5$WEIGHTED_FST, na.rm = TRUE)


# Add category column
fst_highlighted <- fst_windowed %>%
  mutate(FST1 = case_when(
    WEIGHTED_FST == 1 ~ "FST = 1",
    TRUE ~ "Other"
  ))

fst_highlighted <- fst_highlighted %>%
  mutate(Top_1p = case_when(
    WEIGHTED_FST >= fst_threshold_1p ~ "Top 1% FST",
    TRUE ~ "Other"
  ))

fst_highlighted <- fst_highlighted %>%
  mutate(Top_5p = case_when(
    WEIGHTED_FST >= fst_threshold_5p ~ "Top 5% FST",
    TRUE ~ "Other"
  ))

# Plot using manhattan()
# highlight only FST = 1 windows
manhattan(
  fst_highlighted,
  chr = "CHR",
  bp = "BP",
  p = "WEIGHTED_FST",
  snp = "SNP",
  logp = FALSE,
  highlight = fst_highlighted$SNP[fst_highlighted$FST1 == "FST = 1"],
  main = "FST Manhattan Plot with FST = 1 Windows Highlighted",
  ylim = c(0, 1.01)
)

# highlight top 1% FST windows
manhattan(
  fst_highlighted,
  chr = "CHR",
  bp = "BP",
  p = "WEIGHTED_FST",
  snp = "SNP",
  logp = FALSE,
  highlight = fst_highlighted$SNP[fst_highlighted$Top_1p == "Top 1% FST"],
  main = "FST Manhattan Plot with Top 1% FST Windows Highlighted",
  ylim = c(0, 1.01)
)

# highlight top 5% FST windows
manhattan(
  fst_highlighted,
  chr = "CHR",
  bp = "BP",
  p = "WEIGHTED_FST",
  snp = "SNP",
  logp = FALSE,
  highlight = fst_highlighted$SNP[fst_highlighted$Top_5p == "Top 5% FST"],
  main = "FST Manhattan Plot with Top 5% FST Windows Highlighted",
  ylim = c(0, 1.01)
)

# make a table of specific SNPs in the windows with FST = 1
# Initialize list to store results
matched_snps <- list()

# Loop through each window in fst_one
for (i in 1:nrow(fst_one)) {
  chrom <- fst_one$CHROM[i]
  start <- fst_one$BIN_START[i]
  end <- fst_one$BIN_END[i]
  
  # Filter SNPs in this window
  snps_in_window <- fst %>%
    filter(CHROM == chrom & POS >= start & POS <= end)
  
  # Store with window info
  matched_snps[[i]] <- snps_in_window %>%
    mutate(WINDOW_ID = paste0(chrom, ":", start, "-", end))
}


# Combine all results into one data frame
all_matched_snps <- bind_rows(matched_snps)
# save as csv
write_csv(all_matched_snps, "snps_in_fst_one_windows.csv")

# same for top 1%
# Initialize list to store results
matched_snps1 <- list()

# Loop through each window in fst_one
for (i in 1:nrow(fst_top_1)) {
  chrom <- fst_top_1$CHROM[i]
  start <- fst_top_1$BIN_START[i]
  end <- fst_top_1$BIN_END[i]
  
  # Filter SNPs in this window
  snps_in_window <- fst %>%
    filter(CHROM == chrom & POS >= start & POS <= end)
  
  # Store with window info
  matched_snps1[[i]] <- snps_in_window %>%
    mutate(WINDOW_ID = paste0(chrom, ":", start, "-", end))
}


# Combine all results into one data frame
all_matched_snps1 <- bind_rows(matched_snps1)
# save as csv
write_csv(all_matched_snps1, "snps_in_fst_top_1p_windows.csv")


# same for top 5%
# Initialize list to store results
matched_snps5 <- list()

# Loop through each window in fst_one
for (i in 1:nrow(fst_top_5)) {
  chrom <- fst_top_5$CHROM[i]
  start <- fst_top_5$BIN_START[i]
  end <- fst_top_5$BIN_END[i]
  
  # Filter SNPs in this window
  snps_in_window <- fst %>%
    filter(CHROM == chrom & POS >= start & POS <= end)
  
  # Store with window info
  matched_snps5[[i]] <- snps_in_window %>%
    mutate(WINDOW_ID = paste0(chrom, ":", start, "-", end))
}


# Combine all results into one data frame
all_matched_snps5 <- bind_rows(matched_snps5)
# save as csv
write_csv(all_matched_snps5, "snps_in_fst_top_5p_windows.csv")

## plot with ggplot to be able to change colours etsc
# 1) Clean WEIGHTED_FST: set negatives to 0
fst_windowed_clean <- fst_windowed %>%
  mutate(WEIGHTED_FST = pmax(WEIGHTED_FST, 0))

# 2) Recompute threshold after cleaning (top 1%; change to 0.95 for top 5%)
threshold <- quantile(fst_windowed_clean$WEIGHTED_FST, probs = 0.95, na.rm = TRUE)

# 3) Prepare cumulative positions for Manhattan-style x-axis
chr_info <- fst_windowed_clean %>%
  select(CHR, BP) %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(tot = cumsum(chr_len) - chr_len)

fst_plot <- fst_windowed_clean %>%
  inner_join(chr_info %>% select(CHR, tot), by = "CHR") %>%
  mutate(
    BP_cum = BP + tot,
    above_thresh = WEIGHTED_FST >= threshold
  )

axis_df <- fst_plot %>%
  group_by(CHR) %>%
  summarise(center = (min(BP_cum) + max(BP_cum)) / 2, .groups = "drop")

# 4) Plot (no negative y-values; highlighted points on top)
ggplot() +
  geom_point(
    data = filter(fst_plot, !above_thresh),
    aes(x = BP_cum, y = WEIGHTED_FST, color = factor(CHR %% 2)),
    size = 1, alpha = 0.8, shape = 16
  ) +
  scale_color_manual(values = c("0" = "#4D4D4D", "1" = "#B3B3B3")) +
  geom_point(
    data = filter(fst_plot, above_thresh),
    aes(x = BP_cum, y = WEIGHTED_FST),
    shape = 21, fill = "white", color = "black", stroke = 0.5, size = 1
  ) +
  scale_x_continuous(labels = axis_df$CHR, breaks = axis_df$center) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Scaffold",
    y = "Weighted FST"
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none")

