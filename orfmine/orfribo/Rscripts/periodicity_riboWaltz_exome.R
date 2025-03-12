library(riboWaltz)
library(stringr)
library(data.table)
library(ggplot2)

# psite() function from ribowaltz package
# Control of NAs added in case there is not enough material to find best_offset
psite_ribowaltz <- function(data, flanking = 6, start = TRUE, extremity = "auto",
                            plot = FALSE, plot_dir = NULL, plot_format = "png", cl = 99,
                            txt = FALSE, txt_file = NULL) {
  
  if (txt == T | txt == TRUE) {
    options(warn=-1)
    if (length(txt_file) == 0) {
      dir <- getwd()
      txt_file <- paste0(dir, "/duplicates_filtering.txt")
    } else {
      txt_file_split <- strsplit(txt_file, "/")[[1]]
      txt_dir <- paste(txt_file_split[-length(txt_file_split)], collapse = "/")
      if (!dir.exists(txt_dir)) {
        dir.create(txt_dir, recursive = TRUE)
      }
    }
    options(warn=0)
    
    cat("sample\texremity\toffset(nts)\n", file = txt_file)
  }
  
  names <- names(data)
  offset <- NULL
  for (n in names) { 
    cat(sprintf("processing %s\n", n))
    dt <- data[[n]]
    lev <- sort(unique(dt$length))
    if(start == T | start == TRUE){
      base <- 0
      dt[, site_dist_end5 := end5 - cds_start]
      dt[, site_dist_end3 := end3 - cds_start]
    } else {
      base <- -5
      dt[, site_dist_end5 := end5 - cds_stop - base]
      dt[, site_dist_end3 := end3 - cds_stop - base]
    }
    site_sub <- dt[site_dist_end5 <= -flanking & site_dist_end3 >= flanking - 1]
    minlen <- min(site_sub$length)
    maxlen <- max(site_sub$length)
    t <- table(factor(site_sub$length, levels = lev))
    
    # offset
    offset_temp <- data.table(length = as.numeric(as.character(names(t))), percentage = (as.vector(t)/sum(as.vector(t))) * 100)
    offset_temp[, around_site := "T"
    ][percentage == 0, around_site := "F"]
    tempoff <- function(v_dist){
      ttable <- sort(table(v_dist), decreasing = T)
      ttable_sr <- ttable[as.character(as.numeric(names(ttable))+1)]
      ttable_sl <- ttable[as.character(as.numeric(names(ttable))-1)]
      tsel <- rowSums(cbind(ttable > ttable_sr, ttable > ttable_sl), na.rm = T)
      return(as.numeric(names(tsel[tsel == 2][1])))
    }
    
    offset_temp5 <- site_sub[, list(offset_from_5 = tempoff(.SD$site_dist_end5)), by = length]
    offset_temp3 <- site_sub[, list(offset_from_3 = tempoff(.SD$site_dist_end3)), by = length]
    
    # Test to put default values to psite offsets if there is not enough material 
    if(length(unique(offset_temp5$offset_from_5)) < 2) {
      print("There is not enough material for riboWaltz to find the 5prime offset.")
      print("Putting 12nt as 5prime offset to continue the pipeline.")
      offset_temp5$offset_from_5 <- rep.int(12, length(offset_temp5$offset_from_5))
    }
    if(length(unique(offset_temp3$offset_from_3)) < 2) {
      print("There is not enough material for riboWaltz to find the 3prime offset.")
      print("Putting a 3prime offset depending on 5prime offset to continue the pipeline.")
      offset_temp3$offset_from_3 <- offset_temp5$length-offset_temp5$offset_from_5-1
    }
    
    merge_allx <- function(x, y) merge(x, y, all.x = TRUE, by = "length")
    offset_temp  <-  Reduce(merge_allx, list(offset_temp, offset_temp5, offset_temp3))
    
    # adjusted offset
    adj_off <- function(dt_site, dist_site, add, bestoff){
      temp_v <- dt_site[[dist_site]]
      t <- table(factor(temp_v, levels = seq(min(temp_v) - 2, max(temp_v) + add)))
      t[1:2] <- t[3] + 1
      locmax <- as.numeric(as.character(names(t[which(diff(sign(diff(t))) == -2)]))) + 1
      adjoff <- locmax[which.min(abs(locmax - bestoff))]
      ifelse(length(adjoff) != 0, adjoff, bestoff)
    }
    
    best_from5_tab <- offset_temp[!is.na(offset_from_5), list(perc = sum(percentage)), by = offset_from_5
    ][perc == max(perc)]
    best_from3_tab <- offset_temp[!is.na(offset_from_5), list(perc = sum(percentage)), by = offset_from_3
    ][perc == max(perc)]
    
    if(extremity == "auto" &
       ((best_from3_tab[, perc] > best_from5_tab[, perc] &
         as.numeric(best_from3_tab[, offset_from_3]) <= minlen - 2) |
        (best_from3_tab[, perc] <= best_from5_tab[, perc] &
         as.numeric(best_from5_tab[, offset_from_5]) > minlen - 1)) |
       extremity == "3end"){
      best_offset <- as.numeric(best_from3_tab[, offset_from_3])
      line_plot <- "3end"
      adj_tab <- site_sub[, list(corrected_offset_from_3 = adj_off(.SD, "site_dist_end3", 0, best_offset)), by = length]
      offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
      offset_temp[is.na(corrected_offset_from_3), corrected_offset_from_3 := best_offset
      ][, corrected_offset_from_5 := -corrected_offset_from_3 + length - 1]
    } else {
      if(extremity == "auto" &
         ((best_from3_tab[, perc] <= best_from5_tab[, perc] &
           as.numeric(best_from5_tab[, offset_from_5]) <= minlen - 1) |
          (best_from3_tab[, perc] > best_from5_tab[, perc] &
           as.numeric(best_from3_tab[, offset_from_3]) > minlen - 2)) |
         extremity == "5end"){
        best_offset <- as.numeric(best_from5_tab[, offset_from_5])
        line_plot <- "5end"
        adj_tab <- site_sub[, list(corrected_offset_from_5 = adj_off(.SD, "site_dist_end5", 1, best_offset)), by = length]
        offset_temp <- merge(offset_temp, adj_tab, all.x = TRUE, by = "length")
        offset_temp[is.na(corrected_offset_from_5), corrected_offset_from_5 := best_offset
        ][, corrected_offset_from_5 := abs(corrected_offset_from_5)
        ][, corrected_offset_from_3 := abs(corrected_offset_from_5 - length + 1)]
      }
    }
    
    cat(sprintf("best offset: %i nts from the %s\n", abs(best_offset), gsub("end", "' end", line_plot)))
    
    if(txt == T | txt == TRUE){
      cat(sprintf("%s\t%s\t%i\n", n, gsub("end", "'end", line_plot), abs(best_offset)), file = txt_file, append = TRUE)
    }
    
    t <- table(factor(dt$length, levels = lev))
    offset_temp[!is.na(offset_from_5), offset_from_5 := abs(offset_from_5)
    ][, total_percentage := as.numeric(format(round((as.vector(t)/sum(as.vector(t))) * 100, 3), nsmall=4))
    ][, percentage := as.numeric(format(round(percentage, 3), nsmall=4))
    ][, sample := n]
    
    setcolorder(offset_temp, c("length", "total_percentage", "percentage", "around_site", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
    if(start == TRUE | start == T){
      setnames(offset_temp, c("length", "total_percentage", "start_percentage", "around_start", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
      xlab_plot<-"Distance from start (nt)"
    } else {
      setnames(offset_temp, c("length", "total_percentage", "stop_percentage", "around_stop", "offset_from_5", "offset_from_3", "corrected_offset_from_5", "corrected_offset_from_3", "sample"))
      xlab_plot<-"Distance from stop (nt)"
    }
    
    # plot
    if (plot == T | plot == TRUE) {
      options(warn=-1)
      if (length(plot_dir) == 0) {
        dir <- getwd()
        plot_dir <- paste(dir, "/offset_plot", sep = "")
      }
      if (!dir.exists(plot_dir)) {
        dir.create(plot_dir)
      }
      minlen <- ceiling(quantile(site_sub$length, (1 - cl/100)/2))
      maxlen <- ceiling(quantile(site_sub$length, 1 - (1 - cl/100)/2))
      for (len in minlen:maxlen) {
        progress <- ceiling(((len + 1 - minlen)/(maxlen - minlen + 1)) * 25)
        cat(sprintf("\rplotting   %s\r", paste(paste(rep(c(" ", "<<", "-"),
                                                         c(25 - progress, 1, progress)), collapse = ""), " ", as.character(progress*4),
                                               "% ", paste(rep(c("-", ">>", " "), c(progress, 1, 25 - progress)), collapse = ""), sep = "")))
        site_temp <- dt[site_dist_end5 %in% seq(-len + 1, 0) & length == len]
        site_tab5 <- data.table(table(factor(site_temp$site_dist_end5, levels = (-len + 1) : (len))))
        site_temp <- dt[site_dist_end3 %in% seq(0, len - 2) & length == len]
        site_tab3 <- data.table(table(factor(site_temp$site_dist_end3, levels = (-len) : (len - 2))))
        setnames(site_tab5, c("distance", "reads"))
        setnames(site_tab3, c("distance", "reads"))
        site_tab5[, distance := as.numeric(as.character(site_tab5$distance))
        ][, extremity := "5' end"]
        site_tab3[, distance := as.numeric(as.character(site_tab3$distance))
        ][, extremity := "3' end"]
        final_tab <- rbind(site_tab5[distance <= 0], site_tab3[distance >= 0])
        final_tab[, extremity := factor(extremity, levels = c("5' end", "3' end"))]
        
        p <- ggplot(final_tab, aes(distance, reads, color = extremity)) +
          geom_line() +
          geom_vline(xintercept = seq(floor(min(final_tab$distance)/3) * 3, floor(max(final_tab$distance)/3) * 3, 3), linetype = 2, color = "gray90") +
          geom_vline(xintercept = 0, color = "gray50") +
          geom_vline(xintercept = - offset_temp[length == len, offset_from_5], color = "#D55E00", linetype = 2, size = 1.1) +
          geom_vline(xintercept = offset_temp[length == len, offset_from_3], color = "#56B4E9", linetype = 2, size = 1.1) +
          geom_vline(xintercept = - offset_temp[length == len, corrected_offset_from_5], color = "#D55E00", size = 1.1) +
          geom_vline(xintercept = offset_temp[length == len, corrected_offset_from_3], color = "#56B4E9", size = 1.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - len, xmax = -flanking , fill = "#D55E00", alpha = 0.1) +
          annotate("rect", ymin = -Inf, ymax = Inf, xmin = flanking - 1 , xmax = len - flanking - 1, fill = "#56B4E9", alpha = 0.1) +
          labs(x = xlab_plot, y = "Number of read extremities", title = paste(n, " - length=", len, " nts", sep = ""), color= "Extremity") +
          theme_bw(base_size = 20) +
          scale_fill_discrete("") +
          theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), strip.placement = "outside") +
          theme(plot.title = element_text(hjust = 0.5))
        
        if(line_plot == "3end"){
          p <- p + geom_vline(xintercept = best_offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset - len + 1, color = "black", linetype = 3, size = 1.1)
        } else {
          p <- p + geom_vline(xintercept = best_offset, color = "black", linetype = 3, size = 1.1) +
            geom_vline(xintercept = best_offset + len - 1, color = "black", linetype = 3, size = 1.1)
        }
        
        p <- p + 
          scale_x_continuous(limits = c(min(final_tab$distance), max(final_tab$distance)),
                             breaks = seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5), 
                             labels = as.character(seq(floor(min(final_tab$distance)/5) * 5, floor(max(final_tab$distance)/5) * 5, 5) + base))
        
        subplot_dir <- paste(plot_dir, n, sep = "/")
        dir.create(subplot_dir)
        ggsave(paste(subplot_dir, "/", len, ".", plot_format, sep = ""), plot = p, width = 15, height = 5, units = "in")
      }
      cat(sprintf("\rplotting   %s\n",
                  paste(paste(rep(c(" ", "<<", "-"), c(25 - progress, 1, progress)), collapse = ""), " ", 
                        as.character(progress*4), "% ", 
                        paste(rep(c("-", ">>", " "), c(progress, 1, 25 - progress)), collapse = ""), sep = "")))
      options(warn=0)
    }
    
    dt[, c("site_dist_end5", "site_dist_end3") := NULL]
    offset <- rbind(offset, offset_temp)
  }
  return(offset)
}



#result = tryCatch({

#args <- commandArgs(trailingOnly=TRUE)

#config <- args[1]
#gtf_file <- args[2]
#bam_folder <- args[3]
#outdir <- args[4] 

#dir.create(file.path(outdir), showWarnings=F, recursive=TRUE)


# Read config file
#params <- scan(file = file.path(config),
#               what = "character",
#               sep = ":"
#)
result = tryCatch({
  
  args <- commandArgs(trailingOnly=TRUE)
  
#  config <- args[1]
  gtf_file <- args[1]
  bam_folder <- args[2]
  min_size <- as.numeric(args[3])
  max_size <- as.numeric(args[4])
  outdir <- args[5]
  dir.create(file.path(outdir), showWarnings=F, recursive=TRUE)
  
  # Read config file
  #lines <- readLines(config)
  #params <- list()
  #for (line in lines) {
  #  parts <- strsplit(line, ":")[[1]]
  #  if (length(parts) == 2) {
  #    key <- trimws(parts[1])
  #    value <- trimws(parts[2])
  #    params[[key]] <- value
  #  }
  #}

  # Creates annotation table by transcript names
  annotation_db <- riboWaltz::create_annotation(gtf_file)
  annotation_db_transcript_with_cds0l <- data.table(annotation_db)
  annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0,]
  
  # Free unused memory
  rm(list=c("annotation_db","annotation_db_transcript_with_cds0l"))
  gc()

  # Bam files to be computed
  bam_list <- list.files(bam_folder, pattern = "\\.bam$")
  samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
  names(samples) <- str_remove(bam_list, ".bam")
  samples
  
  # Convert BAM to list and filter by fragment length
  reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)
  
  # Filter reads by fragment size
  for (sample in names(reads_list)) {
    reads_list[[sample]] <- reads_list[[sample]][reads_list[[sample]]$length >= min_size & reads_list[[sample]]$length <= max_size, ]
  }
  
  # p-site calculation
  psite_offset <- psite_ribowaltz(reads_list,
                                  flanking = 6,
                                  start = TRUE,
                                  extremity = "auto",
                                  plot = TRUE,
                                  plot_dir = outdir,
                                  plot_format = "tiff",
                                  cl = 100,
                                  txt = TRUE,
                                  txt_file = file.path(outdir, "best_offset.txt")
  )

  reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)

  write.table(psite_offset, file.path(outdir, "psite_offset.csv"), quote = F, row.names = F, sep ="\t")

})

