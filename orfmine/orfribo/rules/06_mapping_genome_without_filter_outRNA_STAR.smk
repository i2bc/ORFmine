rule index_reference_BOWTIE2:
    input: 
         fasta= str(FASTA_PATH) 
    output: 
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    log: 
         str(LOGS_PATH / "Mapping" / "Genome" / "Bowtie2" / "Index" / "Genome_index.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Genome" / "Bowtie2" / "Index" / "Genome_index.benchmark.txt") 
    params: 
         str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Bowtie2" / "Index" / "index_bowtie2")
    threads:
         THREADS_NB
    shell: 
         "bowtie2-build --threads {threads} {input.fasta} {params} &> {log} ;"


rule index_reference_STAR:
    input:
        str(FASTA_PATH)
    output:
        directory(str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Star" / "Index"))
    log:
        str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "Index" / "Genome_index.log")
    benchmark:
        str(BENCHMARKS_PATH / "Mapping" / "Genome" / "Star" / "Index" / "Genome_index.txt")
    shell:
        "mkdir {output} && "
        "../STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} --genomeSAindexNbases 4 > {log} 2>&1 "
 
       
rule Mapping_Genome_STAR_Bowtie2:
    input: 
       fastq = str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")),    
       index_star = str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Star" / "Index"),
       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    output:
       sam_star = str(DATA_PROCESSING_PATH / "Mapping" /  "Genome" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
       met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1")
    params:
       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Bowtie2" / "Index" / "index_bowtie2"),
       prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Star" / "Results" / "{sample}" / "{sample}_"),
       introns = INTRONS_LENGTH
    threads: 
       THREADS_NB
    log:
       final = str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_Genome_Log.final.out"),
       log = str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_Genome_Log.out"),
       sj = str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_Genome_SJ.out.tab"),
       prog = str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_Genome_Log.progess.out"),
       star = str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_star.out"), 
       bowtie2_out = str(LOGS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_bowie2_mapping.out")
    benchmark:
       str(BENCHMARKS_PATH / "Mapping" / "Genome" / "Star" / "{sample}_STAR_Bowtie2_Mapping_Genome.benchmark.txt")
    shell:
        "../STAR-2.7.11b/bin/Linux_x86_64/STAR --readFilesCommand zcat " 
        " --alignIntronMax {params.introns} "
        " --outReadsUnmapped Fastx "
        " --genomeDir {input.index_star}"
        " --runThreadN {threads} " 
        " --readFilesIn {input.fastq} " 
        " --outFileNamePrefix {params.prefix} 2>> {log.star}; " 
        "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -U {output.met1} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}" 
	        
#rule compressed_unmapped_ALL_RNA:
#   input:
#     met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1")
#   output:
#     met1_compressed = str(DATA_PROCESSING_PATH / "Genome_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz")
#   shell:
#     "gzip -c {input.met1} > {output.met1_compressed}"
