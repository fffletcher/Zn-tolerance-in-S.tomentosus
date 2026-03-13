library(tidyverse)

setwd("~/Desktop/Stomentosus_PopGen/PCA/")

# load data from plink
pca <- read_table("pca_zn20_results.eigenvec", col_names = FALSE) #eigenvec
eigenval <- scan("pca_zn20_results.eigenval") #eigenval

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

# add tolerance data

# Read sample lists
tolerant <- read_lines("~/Desktop/Stomentosus_PopGen/top10.txt")
sensitive <- read_lines("~/Desktop/Stomentosus_PopGen/bottom10.txt")

# Create tolerance status data frame
tolerance_df <- tibble(
  ind = c(tolerant, sensitive),
  tolerance = c(rep("Tolerant", length(tolerant)), rep("Sensitive", length(sensitive)))
)

# Join with PCA data
pca2 <- left_join(pca, tolerance_df, by = "ind")

pca2[pca2$ind=="AATTAGATTG-GCATTTAGCG", "tolerance"] <- "Sensitive"


# Plot PCA with color
ggplot(pca2, aes(PC1, PC2, color = tolerance)) +
  geom_point(size = 3) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_manual(values = c("Tolerant" = "#0072B2", "Sensitive" = "#E69F00"))

# label 220 and 230 on PCA

# Define the isolates to highlight
highlight_df <- tibble(
  ind = c("ACCCACAACT-CGGTGGTAAG", "GTATGGCTTC-CCGATCGCCA"),
  label = c("ZnS", "ZnT")
)
          
# Join with PCA data to get coordinates
highlight_coords <- left_join(highlight_df, pca2, by = "ind")

# Define offset for arrows (adjust as needed)
highlight_coords <- highlight_coords %>%
  mutate(xend = PC1 + 0.25, yend = PC2 + 0.1)

# Plot PCA with arrows and labels
ggplot(pca2, aes(PC1, PC2, color = tolerance)) +
  geom_point(size = 5) +
  geom_segment(data = highlight_coords, aes(x = xend, y = yend, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.1, "cm")), color = "black", alpha = 0.2) +
  geom_text(data = highlight_coords, aes(x = xend, y = yend, label = label),
            hjust = -0.1, vjust = .5, size = 3,  alpha = 1, show.legend = F) +
  coord_equal() +
  theme_light() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_manual(values = c("Tolerant" = "#0072B2", "Sensitive" = "#E69F00")) +
  labs(color = "Tolerance Status")

# make 220 and 230 open ciricles instead
# Define isolates to highlight (open circles)
highlight_ids <- c("ACCCACAACT-CGGTGGTAAG", "GTATGGCTTC-CCGATCGCCA")

# Add a new column to indicate shape
pca2 <- pca2 %>%
  mutate(shape_flag = if_else(ind %in% highlight_ids, "highlight", "normal"))

# Plot PCA with custom shapes

ggplot(pca2, aes(PC1, PC2, color = tolerance)) +
  geom_point(aes(shape = shape_flag), size = 5, stroke = 2) + 
  scale_shape_manual(values = c("normal" = 16, "highlight" = 21)) +
  coord_equal() +
  theme_bw() + 
  theme(axis.text  = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_color_manual(values = c("Tolerant" = "#0072B2", "Sensitive" = "#E69F00")) +
  labs(color = "Zn Tolerance", shape = NULL) +
  geom_segment(data = highlight_coords, aes(x = xend, y = yend, xend = PC1, yend = PC2), 
               color = c("#E69F00", "#0072B2"), alpha = .5, linewidth = 1) +
  geom_text(data = highlight_coords, aes(x = xend, y = yend, label = label),
            hjust = -0.1, vjust = .5, size = 5,  alpha = 1, show.legend = F) +
  guides(shape = "none")

# zoom on this as you can't see 220 unless you do ^^
