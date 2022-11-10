library(riboWaltz)
library(stringr)
library(data.table)


result = tryCatch({

setwd("/workdir/")
  
# Path to working directory
local_path <- "/workdir/orfribo/"
dir.create(paste0(local_path, "RESULTS/ORFribo/riboWaltz/"), showWarnings = F)

# Read config file
params <- scan(file = paste0("/workdir/config.yaml"),
               what = "character",
               sep = ":"
)
args <- commandArgs(trailingOnly=TRUE)

gtf_file <- args[1]

# Creates annotation table by transcript names
annotation_db <- riboWaltz::create_annotation(gtf_file)
annotation_db_transcript_with_cds0l <- data.table(annotation_db)
annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0,]
# Free unused memory
rm(list=c("annotation_db","annotation_db_transcript_with_cds0l"))
gc()

# Bam files to be computed
bam_folder <- paste0(local_path, "RESULTS/ORFribo/BAM_exome/")

bam_list <- list.files(bam_folder, pattern = "\\.bam$")

samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
names(samples) <- str_remove(bam_list, ".bam")
samples

reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)


# p-site calculation
source(paste0("/ORFmine/orfribo/RiboDoc_BAM2Reads/tools/ribowaltz_psite_with_NA_control.R"))
psite_offset <- psite_ribowaltz(reads_list,
                                flanking = 6,
                                start = TRUE,
                                extremity = "auto",
                                plot = TRUE,
                                plot_dir = paste0(local_path, "RESULTS/ORFribo/riboWaltz/"),
                                plot_format = "tiff",
                                cl = 100,
                                txt = TRUE,
                                txt_file = paste0(local_path, "RESULTS/ORFribo/riboWaltz/best_offset.txt")
)


reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)

write.table(psite_offset, paste0(local_path, "RESULTS/ORFribo/riboWaltz/psite_offset.csv"), quote = F, row.names = F, sep ="\t")

})
