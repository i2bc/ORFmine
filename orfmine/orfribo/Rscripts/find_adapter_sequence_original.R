library(stringr)

# Specify the path to the working directory
args <- commandArgs(trailingOnly=TRUE)
local_path <- args[2]

# For file size loaded (entire fastq usually too big for R)
ibs <- 8192
ibs_count <- 50

# Default sequence size checked
Adapt_size <- 10

# Fraction (percentage) of lines which should contain the sequence to be a potential adapter
fraction_threshold <- 0.7

# Take the path and fastq name as argument
fastq <- args[1]
sample <- sub(".fastq.gz", "", basename(fastq))

# Open a part of the fastq file to try finding the adapter sequence corresponding to this sequening
try(File <- system(paste0("gzip -cd ", fastq," 2> /dev/null | dd ibs=",ibs," count=", ibs_count, " | grep -P '^[ATGC]+$'"), intern = TRUE))
File <- as.data.frame(File)

# Number of lines to check (arbitraty)
lines_check=1000
if(nrow(File) < lines_check) {lines_check <- nrow(File)}

# Initialize list of potential adapters
List_adapters <- NULL

# Find potential adapters sequence
for(Ligne in 1 : lines_check)
{
  #Scan line
  for(Position1 in 1 : ((nchar(File[Ligne,1]) - Adapt_size )))
  {
    # Sequence checked in whole file (found at line Ligne from position Position1 with Adapt_size as length)
    Pattern <- str_sub(File[Ligne,1], Position1, (Position1 + Adapt_size))
    
    # Number of this sequence occurrences in all fastq lines selected
    Occurence <- suppressWarnings(str_count(File, Pattern))
    
    # If this sequence is found in less than 70% (0.7 : arbitrary) of the fastq lines : the window is moved
    if(Occurence < (nrow(File) * fraction_threshold)) {next}
    
    # Enlarging the window of the checked sequence
    for(Position2 in (Position1 + Adapt_size) : nchar(File[Ligne,1]))
    {
      # New checked sequence
      Pattern2 <- str_sub(File[Ligne,1], Position1, Position2)
      New_Occurence <- suppressWarnings(str_count(File, Pattern2))
      
      if(New_Occurence < (nrow(File) * fraction_threshold))
      {
        Pattern2 <- str_sub(File[Ligne,1], Position1, Position2-1)
        break
      }
      # If the checked sequence is long enough (more than 20nt (arbitrary))
      if(nchar(Pattern2) == 20) {break}
    }
    
    # Sequence is saved as a potential adapter
    List_adapters <- rbind(List_adapters, Pattern2)
    break
  }
  df <- as.data.frame(table(List_adapters))
  
  if(nrow(df) > 0) {
    df <- df[order(df$Freq, decreasing = T),]
    
    # If more than half the checked lines (arbitrary) have been used to make patterns
    # and if the best checked sequence appears in more than 50% (0.5 : arbitrary) of the fastq lines
    # The checked sequence is kept as a potential adapter
    if(Ligne > (lines_check * 0.5) & (max(df$Freq)/Ligne) > fraction_threshold)
    {
      print(paste0("Adapter sequence found for fastq '", basename(fastq), "' : ", Pattern2))
      write(Pattern2, paste0(local_path, sample, ".txt"), append = TRUE)
      break
    }
  }
}
