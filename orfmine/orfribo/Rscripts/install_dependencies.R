list.of.packages <- c("stringr", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)


# sudo apt install libcurl4-openssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 
# R -e "library('devtools');devtools::install_github('LabTranslationalArchitectomics/riboWaltz@v1.2.0', dependencies = FALSE);"
