# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Update 'Matrix' package
install.packages("Matrix")

# Package names
packages <- c("ggplot2", "ggbeeswarm", "ggpubr", "plotrix", "RColorBrewer", "ineq", "ggExtra", "ggridges", "proxy", "dplyr", "seqinr", "stringr")

# Install packages not yet installed with dependencies
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], dependencies = TRUE)
}

library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(plotrix)
library(RColorBrewer)
library(ineq)
library(ggExtra)
library(ggridges)
library(proxy)
library(dplyr)
library(seqinr)
library(stringr)

coding.ribo.reads <- read.csv(paste("/data/work/I2BC/fadwa.elkhaddar/BIM/ORFMINE/ORFmine/orfmine/orfribo/Extra_Scripts/Coding_ORFs.csv", sep = ""))

reads.cutoff <- 50
cov.cutoff <- 0.3

map.cases.on.triangle <- function(landskape.data, landskape.palette, landskape.alpha = 0.7, mapping.data, by = "black", title, alpha, legend.title, legend.angle = 0, legend.hjust = 0.5, legend.vjust = 0, low = "blue", mid = "grey", high = "red") {
  # ==================================================== #
  # This function generates the triangle. 
  # landskape.data -> Is the dataset of reference on which you want to map specific properites
  # landskape.palette -> The color palette to use for the density of the landskape 
  # mapping.data -> The table of the data to map on top of the landskape dataset with POINTS
  # by -> The property by which to color the points (by default "Black" if you want only to see the localisation of the points)
  # ==================================================== #
  midpoint <- (max(as.numeric(by), na.rm = TRUE) - abs(min(as.numeric(by), na.rm = TRUE))) / 2
  ggplot(data = landskape.data) +  
    stat_density_2d_filled(aes(x = as.numeric(Perc_p0), y = as.numeric(Perc_p1), fill = ..level..), geom = "polygon", alpha = landskape.alpha) + 
    scale_fill_brewer(palette = landskape.palette) +
    geom_polygon(data = data.frame(x = c(0, 100, 0), y = c(0, 0, 100)), aes(x, y), alpha = 0, color = 'black', size = 1.5) +
    geom_point(data = mapping.data, aes(x = as.numeric(Perc_p0), y = as.numeric(Perc_p1), color = by), size = 0.3, alpha = alpha) +
    geom_polygon(data = data.frame(x = c(0, 34, 0), y = c(66, 66, 100)), aes(x, y), alpha = 0, color = 'black', size = 1) +
    geom_polygon(data = data.frame(x = c(66, 100, 66), y = c(0, 0, 34)), aes(x, y), alpha = 0, color = 'black', size = 1) +
    geom_rect(mapping = aes(xmin = 0, xmax = 33, ymin = 0, ymax = 33), color = "black", alpha = 0, size = 1) +
    geom_rect(mapping = aes(xmin = 0, xmax = 33, ymin = 33, ymax = 66), color = "black", alpha = 0, size = 1) +
    scale_color_gradient2(midpoint = midpoint, low = low, mid = mid, high = high, space = "Lab") +
    theme_minimal() + 
    xlab("Fraction of F0 reads") + 
    ylab("Fraction of F1 reads") +
    guides(color = guide_colorbar(title.position = "top", title.hjust = 1, title = legend.title), fill = FALSE) + 
    ggtitle(title) +
    scale_y_continuous(breaks = c(0, 33, 66, 100)) +
    scale_x_continuous(breaks = c(0, 33, 66, 100)) +
    theme(
      legend.title = element_text(size = 7), 
      legend.text = element_text(size = 6, angle = legend.angle, vjust = legend.vjust, hjust = legend.hjust),
      legend.margin = margin(r = 0, l = 0, t = -3, b = -3),
      legend.key.size = unit(0.4, "cm"),
      legend.title.align = 0.5,
      legend.position = c(0.75, 0.7),
      legend.direction = "horizontal",
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title = element_text(size = 9, face = "plain"),
      axis.text.x = element_text(size = 9, angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

##########################
# CDS protein abundance  #
##########################
tab.tmp <- coding.ribo.reads[coding.ribo.reads$Num_reads >= 50 & coding.ribo.reads$coverage >= cov.cutoff & !is.na(coding.ribo.reads$HCA), ]
tab.tmp <- tab.tmp[order(tab.tmp$HCA), ]
coding.triangle.abundance <- map.cases.on.triangle(
  landskape.data = coding.ribo.reads[coding.ribo.reads$Num_reads >= 50 & !is.na(coding.ribo.reads$HCA), ],
  landskape.palette = "Purples", landskape.alpha = 0,
  mapping.data = tab.tmp,
  by = tab.tmp$HCA,
  title = "Abundance", alpha = 1, legend.title = "log10(ppm)",
  low = brewer.pal(n = 9, name = "Oranges")[1], mid = brewer.pal(n = 9, name = "Oranges")[4], high = brewer.pal(n = 9, name = "Oranges")[7]
)
tab.tmp <- NULL

pdf("/data/work/I2BC/fadwa.elkhaddar/BIM/ORFMINE/ORFmine/orfmine/orfribo/Extra_Scripts/FigureS9.pdf", width = 8, height = 4)
print(coding.triangle.abundance)
dev.off()

