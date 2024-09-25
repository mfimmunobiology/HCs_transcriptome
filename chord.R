library(networkD3)
library(circlize)
library(data.table)
library(randomcoloR)
library(dplyr)
library(ComplexHeatmap)

workingdir <- "/home/vivischuch/Documents/chord/"
setwd(workingdir)

# Read the data
df1 <- fread("hallmark_all_edges.csv", header="auto")

# Extract unique values from the 'from' and 'to' columns
from_values <- unique(df1$fom)
to_values <- unique(df1$to)

df1$fom <- factor(df1$fom, levels = c("HIV_EG", "CMV_EG", "HIV_T", "CMV_T"))

# Calculate the degree of each 'to' sector
to_degrees <- df1[, .(degree = sum(abs(value))), by = to]

# Order 'to' values by their degree
# to_values_ordered <- to_degrees[order(-degree)]$to  # Descending order
to_values_ordered <- c(
  "Interferon gamma response", "Interferon alpha response", "Allograft rejection", "Inflammatory response",  # Immune Response
  "IL6 JAK STAT3 signaling", "TNFa signaling via NFKB", "KRAS signaling DN", "KRAS signaling UP",            # Signaling Pathways
  "Apoptosis", "Epithelial mesenchymal transition",                                                          # Apoptosis and Cell Death
  "Reactive oxygen species pathway", "Hypoxia", "Xenobiotic metabolism",                                     # Metabolic and Stress Responses
  "Myogenesis", "Angiogenesis", "Coagulation", "Apical junction",                                            # Structural and Developmental
  "Estrogen response EARLY", "Estrogen response LATE",                                                       # Hormonal Responses
  "UV response DN", "WNT beta catenin signaling"                                                             # Others
)

# Combine all unique sectors into a single order vector, ensuring 'from' are on the left and ordered 'to' on the right
# sector_order <- c(from_values, to_values_ordered)
sector_order <- c(levels(df1$fom), to_values_ordered)


# Define the color palette with the correct number of colors
n <- length(sector_order)
grid.col <- distinctColorPalette(n)
names(grid.col) <- sector_order

# Define color scale based on NES values
# Define the color function with five colors
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1, 2, 3),
  colors = c("blue", "lightblue", "white", "salmon", "red")
)

svg(paste0(workingdir, "chord_diagram.svg"), width = 10, height = 7)  # Define size of the SVG
png(paste0(workingdir, "chord_diagram.png"), width = 1200, height = 800, res = 150)

# Plot the chord diagram with the NES-based color
par(mfrow = c(1, 1))

# Initialize the circos plot with proper sector order
# Initialize the circos plot with proper sector order
circos.par(
  start.degree = -103,  # Rotate the starting point
  gap.after = c(rep(2, length(levels(df1$fom)) - 1), 10, rep(2, length(to_values_ordered) - 1), 10)
)

# Create the chord diagram
chordDiagram(
  df1[, c("fom", "to", "value")],
  grid.col = grid.col,
  link.visible = df1$value > 0 | df1$value < 0,  # Show both positive and negative links
  col = col_fun(df1$value),  # Map NES values to colors
  transparency = 0.5,  # Optional transparency for aesthetic purposes
  annotationTrack = c("names", "grid"),
  annotationTrackHeight = c(0.03, 0.03),
  preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df1))))),
  order = sector_order,  # Use the defined order with ordered 'to' sectors
  directional = 1,       # Set directional to indicate connections from 'from' to 'to'
  direction.type = c("arrows"), # Optional to show arrows
  link.arr.type = "big.arrow"  # Optional to show big arrows
)

# Customize sector labels on the first track
circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "outside",  # Change to "inside" or "reverse.clockwise" for different alignment
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.7
    )
  },
  bg.border = NA
)
# Define the legend with only whole numbers as labels
lgd <- Legend(
  title = "NES Value",
  col_fun = col_fun,
  at = c(-1, 0, 1, 2, 3),  # Set the breaks to whole numbers
  labels = c("-1", "0", "1", "2", "3"),  # Labels corresponding to the breaks
  direction = "vertical"  # Set the direction of the legend
)

# Draw the legend on the left side
draw(lgd, x = unit(2, "cm"), y = unit(4, "cm"), just = c("left", "bottom"))

# Clear circos plot
circos.clear()

dev.off()
