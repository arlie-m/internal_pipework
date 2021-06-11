# For the full dendrogram
library(plyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)

set.seed(10)

# The source data
mat <- matrix(rnorm(24 * 10, mean = 1, sd = 2), 
              nrow = 24, ncol = 10, 
              dimnames = list(paste("g", 1:24, sep = ""), 
                              paste("sample", 1:10, sep = "")))

sample_names <- colnames(mat)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(mat)))
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = sample_names) %>%
  dplyr::mutate(x_center = (1:n()), 
         width = 1)

# Neglecting the gap parameters
heatmap_data <- mat %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

# Heatmap plot
plt_hmap <- ggplot(heatmap_data, 
                   aes(x = x_center, y = y_center, fill = expr, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_fill_gradient2("expr", high = "darkred", low = "darkblue") +
  scale_x_continuous(breaks = sample_pos_table$x_center, 
                     labels = sample_pos_table$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Sample", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank())

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse(expand = c(0, 0.5)) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

library(cowplot)
plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(1, 1))
