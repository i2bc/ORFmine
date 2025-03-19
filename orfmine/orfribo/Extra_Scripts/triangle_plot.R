# Fonction pour vérifier et installer les packages nécessaires
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        install.packages(pkg, repos = "http://cran.rstudio.com/", dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        message(paste("Erreur lors de l'installation de", pkg, ":", e$message))
      })
    }
  }
}

# Mettre à jour le package Matrix si nécessaire
if (packageVersion("Matrix") < "1.6.0") {
  install.packages("Matrix", repos = "http://cran.rstudio.com/")
  library(Matrix)
}

# Installer les dépendances critiques
critical_dependencies <- c("pbkrtest", "MatrixModels", "quantreg", "car", "rstatix")
install_if_missing(critical_dependencies)

# Liste des packages nécessaires
required_packages <- c("ggplot2", "ggbeeswarm", "ggpubr", "plotrix", "RColorBrewer", "ineq", "ggExtra", "ggridges", "proxy", "dplyr", "seqinr", "stringr")
install_if_missing(required_packages)

# Récupération des arguments de la ligne de commande
args <- commandArgs(trailingOnly = TRUE)

# Affectation des arguments
input_file <- args[1]       # Fichier CSV d'entrée
output_file <- args[2]      # Fichier PDF de sortie
feature <- args[3]          # Nom de la colonne à utiliser pour la colorisation
low_color <- ifelse(length(args) >= 4, args[4], "blue")   # Couleur basse par défaut "blue"
mid_color <- ifelse(length(args) >= 5, args[5], "grey")   # Couleur moyenne par défaut "grey"
high_color <- ifelse(length(args) >= 6, args[6], "red")   # Couleur haute par défaut "red"

# Paramètres de filtrage
reads.cutoff <- 50
cov.cutoff <- 0.3

# Fonction pour générer le triangle
map.cases.on.triangle <- function(landskape.data, landskape.alpha = 0.7, mapping.data, by, title, alpha, legend.title, low = "blue", mid = "grey", high = "red") {
  midpoint <- (max(as.numeric(by), na.rm = TRUE) - abs(min(as.numeric(by), na.rm = TRUE))) / 2
  ggplot(data = landskape.data) +  
    stat_density_2d_filled(aes(x = as.numeric(Perc_p0), y = as.numeric(Perc_p1), fill = ..level..), geom = "polygon", alpha = landskape.alpha) + 
    geom_polygon(data = data.frame(x = c(0, 100, 0), y = c(0, 0, 100)), aes(x, y), alpha = 0, color = 'black', size = 1.5) +
    geom_point(data = mapping.data, aes(x = as.numeric(Perc_p0), y = as.numeric(Perc_p1), color = by), size = 0.3, alpha = alpha) +
    scale_color_gradient2(midpoint = midpoint, low = low, mid = mid, high = high, space = "Lab") +
    theme_minimal() + 
    xlab("Fraction of F0 reads") + 
    ylab("Fraction of F1 reads") +
    guides(color = guide_colorbar(title.position = "top", title.hjust = 1, title = legend.title), fill = FALSE) + 
    ggtitle(title)
}

# Lecture des données
coding.ribo.reads <- read.csv(input_file)

# Vérifier si la colonne feature existe
if (!(feature %in% colnames(coding.ribo.reads))) {
  stop(paste("Erreur : la colonne", feature, "n'existe pas dans le fichier d'entrée."))
}

# Filtrage des données
tab.tmp <- coding.ribo.reads[coding.ribo.reads$Num_reads >= reads.cutoff & coding.ribo.reads$coverage >= cov.cutoff & !is.na(coding.ribo.reads[[feature]]), ]
tab.tmp <- tab.tmp[order(tab.tmp[[feature]]), ]

# Génération du graphique
coding.triangle.abundance <- map.cases.on.triangle(
  landskape.data = coding.ribo.reads[coding.ribo.reads$Num_reads >= reads.cutoff & !is.na(coding.ribo.reads[[feature]]), ],
  landskape.alpha = 0,
  mapping.data = tab.tmp,
  by = tab.tmp[[feature]],
  title = feature, alpha = 1, legend.title = "log10(ppm)",
  low = low_color,
  mid = mid_color,
  high = high_color
)

# Sauvegarde du graphique
pdf(output_file, width = 8, height = 4)
print(coding.triangle.abundance)
dev.off()
