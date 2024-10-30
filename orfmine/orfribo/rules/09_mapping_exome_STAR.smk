rule index_Exome_BOWTIE2:
    input: 
         fasta = str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.nfasta")
    output: 
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    log: 
         str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "Exome_index.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "Exome_index.benchmark.txt")
    params: 
         str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2")
    threads:
         THREADS_NB
    shell: 
         "bowtie2-build --threads {threads} {input.fasta} {params} &> {log} ;"
	 
	 
rule index_Exome_STAR:
    input:
        str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.nfasta")
    output:
        directory(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Index"))
    log:
        str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Index" / "Exome_index.log")
    benchmark:
        str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Star" / "Index" / "Exome_index.txt")
    shell:
        "mkdir {output} && "
        "STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} --genomeSAindexNbases 7 > {log} 2>&1"
 
if config.get('rna_to_exclude'):       
	rule Mapping_Exome_STAR_Bowtie2:
	    input: 
	       fastq= str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}" / "{sample}_Unmapped.fastq.gz"),    
	       index_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Index"),
	       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
	    output:
	       sam_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
	       met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1"),
	       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
               log = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Log.final.out")
	    params:
	       prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_"),
	       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2"),
               multi_map = MULTIMAPPING,
               introns = INTRONS_LENGTH
	    log:
	       final = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_Exome_Log.final.out"),
	       log = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_Exome_Log.out"),
	       sj = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_Exome_SJ.out.tab"),
	       prog = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_Exome_Log.progess.out"),
	       star = str(LOGS_PATH / "Mapping" / "Exome" / "Star" /"Results" / "{sample}" / "Exome_{sample}_star.out"), 
	       bowtie2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Exome_{sample}_bowtie2_star_mapping.txt")
	    benchmark:
	       str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_STAR_Bowtie2_Mapping_Exome.benchmark.txt")
	    shell:
	       "STAR --readFilesCommand zcat " 
	       " --alignIntronMax {params.introns} "
	       " --outReadsUnmapped Fastx "
	       " --genomeDir {input.index_star}"
	       " --runThreadN 20 " 
	       " --readFilesIn {input.fastq} " 
	       " --outFileNamePrefix {params.prefix}" 
	       " --outFilterMultimapNmax {params.multi_map} > {log.star} 2>&1 ; "
	       "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -k {params.multi_map} -U {output.met1} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}"

else: 
	rule Mapping_Exome_STAR_Bowtie2:
	    input: 
	       fastq= str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")),    
	       index_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Index"),
	       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
	    output:
               sam_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
               met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1"),
               sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
               log = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Log.final.out")
	    params:
               prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_"),
               index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2"),
               multi_map = MULTIMAPPING,
               introns = INTRONS_LENGTH
	    log:
               final = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.final.out"),
               log = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.out"),
               sj = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_SJ.out.tab"),
               prog = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.progess.out"),
               star = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Exome_{sample}_star.out"),
               bowtie2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Exome_{sample}_bowtie2_star_mapping.txt")
	    benchmark:
               str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_STAR_Bowtie2_Mapping_Exome.benchmark.txt")
	    shell:
               "STAR --readFilesCommand zcat "
               " --alignIntronMax {params.introns} "
               " --outReadsUnmapped Fastx "
               " --genomeDir {input.index_star}"
               " --runThreadN 20 "
               " --readFilesIn {input.fastq} "
               " --outFileNamePrefix {params.prefix}"
               " --outFilterMultimapNmax {params.multi_map} > {log.star} 2>&1 ; "
               "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -k {params.multi_map} -U {output.met1} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}"

#rule compressed_unmapped_Exome:
#    input: 
#      met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1")
#    output: 
#      met1_compressed = str(DATA_PROCESSING_PATH / "Exome_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz")
#    shell:
#      "gzip -c {input.met1} > {output.met1_compressed}"
      
