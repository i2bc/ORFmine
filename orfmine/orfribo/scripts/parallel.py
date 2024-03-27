import BAM2Reads
import os 



gff = "/home/fadwa.el-khaddar/bam2reads/mapping_orf_Scer.gff"

def gff_chunk(gff):
	with open(gff, "r") as f:
		temp_file= None
		for line in f:
			if lines_written % chunk_size == 0:
				if temp_file:
					temp_file.close()
				temp_file = open(f"temp_file_{file_count}.gff", "w")
				file_count += 1
			temp_file.write(line)
			lines_written += 1
		if temp_file:
			temp_file.close()
		 
