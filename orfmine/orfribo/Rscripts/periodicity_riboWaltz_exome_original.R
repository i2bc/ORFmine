library(riboWaltz)
library(stringr)
library(data.table)


result = tryCatch({
  
  args <- commandArgs(trailingOnly=TRUE)
  
  config <- args[1]
  gtf_file <- args[2]
  bam_folder <- args[3]
  outdir <- args[4]
  
  dir.create(file.path(outdir), showWarnings=F, recursive=TRUE)
  
  # Read config file
  lines <- readLines(config)
  params <- list()
  for (line in lines) {
    parts <- strsplit(line, ":")[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      value <- trimws(parts[2])
      params[[key]] <- value
    }
  }
  


# Creates annotation table by transcript names
annotation_db <- riboWaltz::create_annotation(gtf_file)
annotation_db_transcript_with_cds0l <- data.table(annotation_db)
annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0,]
# Free unused memory
rm(list=c("annotation_db","annotation_db_transcript_with_cds0l"))
gc()

# Bam files to be computed
#bam_folder <- paste0(local_path, "ORFribo/DATA_PROCESSING/SAMtools/Exome")

bam_list <- list.files(bam_folder, pattern = "\\.bam$")

samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
names(samples) <- str_remove(bam_list, ".bam")
samples

reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)


# p-site calculation
source(paste0("/data/work/I2BC/fadwa.elkhaddar/BIM/ORFMINE/ORFmine/orfmine/orfribo/Rscripts/ribowaltz_psite_with_NA_control.R"))
psite_offset <- psite_ribowaltz(reads_list,
                                flanking = 6,
                                start = TRUE,
                                extremity = "auto",
                                plot = TRUE,
                                plot_dir = paste0(outdir),
                                plot_format = "tiff",
                                cl = 100,
                                txt = TRUE,
                                txt_file = paste0(outdir, "best_offset.txt")
)


reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)

write.table(psite_offset, paste0(outdir, "psite_offset.csv"), quote = F, row.names = F, sep ="\t")


})
