library(VennDiagram)
library(data.table)
library(grid)
library(scales)
setwd("~/Desktop/compGenome2/")

pres <- fread("results/core_accessory_singleton/orthogroup_presence.tsv")

if (ncol(pres) < 3) {
  stop("orthogroup_presence.tsv must have at least 3 columns: Orthogroup and 2 isolates")
}

orth_col <- colnames(pres)[1]
iso_cols  <- colnames(pres)[2:3]

label_map <- setNames(c("ZnT", "ZnS"), iso_cols)
friendly  <- unname(label_map[iso_cols])

iso_colors <- c("ZnS" = "#E69F00",
                "ZnT" = "#0072B2")
fill_cols  <- unname(iso_colors[friendly])

set1 <- pres[get(iso_cols[1]) == 1, get(orth_col)]
set2 <- pres[get(iso_cols[2]) == 1, get(orth_col)]

n1  <- length(set1)
n2  <- length(set2)
n12 <- length(intersect(set1, set2))

if (!dir.exists("results")) dir.create("results", recursive = TRUE)

venn_file <- file.path("results", "venn_orthogroups.png")
png(venn_file, width = 1800, height = 1400, res = 200)
grid.newpage()

venn <- draw.pairwise.venn(
  area1      = n1,
  area2      = n2,
  cross.area = n12,
  scaled     = FALSE,
  category   = friendly,
  fill       = fill_cols,
  alpha      = c(0.5, 0.5),
  lty        = "solid",
  lwd        = 2,
  cex        = 1.2,
  cat.cex    = 1.4,
  cat.col    = fill_cols,
  cat.pos    = c(-45, 45),
  cat.dist   = c(0.03, 0.03)
)

grid.draw(venn)
dev.off()
message(sprintf("Wrote Venn: %s", venn_file))
