# For the full dendrogram
library(plyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(cowplot)
library(wesanderson)
library(RColorBrewer)

species_names <- colnames(community_matrix)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(vegdist(community_matrix, method = "bray"))) %>% 
  set("leaves_pch", 19) %>% set("leaves_col", sample_pos_table$labs_colour)
dend_data <- dendro_data(dend)

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))
# Extract dendrogram terminal segment data
dendrogram_ends <- with(segment(dend_data) %>%
  filter(yend == 0),
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the sample labels
sample_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, id = as.character(label), height = 1)) %>% 
  left_join(site_info, by = "id")
# bind the sample information to the dendrogram ends for plotting
dendrogram_ends <- bind_cols(dendrogram_ends, sample_pos_table)

# creating custom colours for the labels of stages
axiscolour = brewer_pal(type = "seq", palette = "Dark2")(6)
#show_col(brewer_pal(type = "seq", palette = "Dark2")(6))
#show_col(viridis_pal(direction = -1, begin = 0.2, end = 0.85)(4))
#viridis_labs = viridis_pal(direction = -1, begin = 0.2, end = 0.85)(4)
sample_pos_table <- sample_pos_table %>% 
  mutate(labs_colour = case_when(stage == 1 ~ axiscolour[1],
                                 stage == 2 ~ axiscolour[2],
                                 stage == 3 ~ axiscolour[3],
                                 stage == 4 ~ axiscolour[6]))
dend_palette <- axiscolour[c(1, 2, 3, 6)]
# Table to position the samples
species_pos_table <- data.frame(taxon = species_names) %>%
  dplyr::mutate(x_center = (1:n()), 
         width = 1)

# Neglecting the gap parameters
heatmap_data <- community_matrix %>% 
  mutate("sample_id" = rownames(community_matrix)) %>% 
  pivot_longer(cols = !sample_id, names_to = "taxon", values_to = "density")%>%
  left_join(sample_pos_table) %>%
  left_join(species_pos_table)

# changing 0 to NA so that samples/samples with low density appear different to those with none.
heatmap_data[heatmap_data == 0] <- NA

# Limits for the vertical axes
sample_axis_limits <- with(
  sample_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

# Heatmap plot
plt_hmap <- ggplot(heatmap_data, 
                   aes(x = x_center, y = y_center, fill = density, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_fill_viridis_c("Density", option = "cividis", na.value = "white", begin = 0.2) +
  scale_x_continuous(breaks = species_pos_table$x_center, 
                     labels = species_pos_table$taxon,
                     expand = c(0, 0)
                     ) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = sample_pos_table[, "y_center"], 
                     labels = rep("", nrow(sample_pos_table)),
                     limits = sample_axis_limits,
                     expand = c(0, 0)
                     ) + 
  labs(x = "Taxa", y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  geom_point(data = dendrogram_ends, aes(x= xend, y = yend, colour = as_factor(sample_pos_table$stage))) +
  scale_x_reverse() + 
  scale_y_continuous(breaks = sample_pos_table$y_center, 
                     labels = sample_pos_table$id, 
                     limits = sample_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "Sample", colour = sample_pos_table$stage, size = "") +
  scale_colour_manual("Stage", values = dend_palette) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(color=sample_pos_table$labs_colour, face = "bold"),
        legend.position = "left",
        legend.direction = "vertical"
        
        )

plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(1, 1.5))
