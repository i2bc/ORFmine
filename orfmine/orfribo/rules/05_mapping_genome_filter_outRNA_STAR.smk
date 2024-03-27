rule index_reference_BOWTIE2:
    input: 
         fasta= str(FASTA_PATH)
    output: 
         expand(str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2.{extb}.bt2",extb=BOWTIE2))
    log: 
         str(LOGS_PATH / "Mapping" / "Orfium"/ "Bowtie2" / "Index"/ "Orfium_index.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Orfium"/ "Bowtie2" / "Index" / "Orfium_index.benchmark.txt") 
    params: 
         str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2")
    threads:
         THREADS_NB
    shell: 
         "bowtie2-build --threads {threads} {input.fasta} {params} &> {log} ;"


rule index_reference_STAR:
    input:
        str(FASTA_PATH)
    #params:
    #    conda_env = "star"
    output:
        directory(str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Star" / "Index"))
    log:
        str(LOGS_PATH / "Mapping" / "Orfium" / "Star" / "Index" / "Orfium_index.log")
    benchmark:
        str(BENCHMARKS_PATH / "Mapping" / "Orfium" / "Star" / "Index" / "Orfium_index.txt")
    shell:
        "set +eu && "
        ". $(conda info --base)/etc/profile.d/conda.sh && "
        "conda activate {params.conda_env} && "
        "mkdir {output} && "
        "STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} --genomeSAindexNbases 4"
 
if config.get('fasta_outRNA'):       
	rule Mapping_Genome_STAR_Bowtie2:
	    input: 
	       fastq = str(DATA_PROCESSING_PATH / "Trimmed_Filtred_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz"),    
	       index_star = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Star" / "Index"),
	       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2.{extb}.bt2",extb=BOWTIE2))
	    output:
	       sam_star = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
	       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
	       met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1")
	    params:
	       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2"),
	       prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Star" / "Results"/ "{sample}" / "{sample}_"),
	       conda_env = "star"
	    threads: 
	       config["multi_threads_nbr"]
	    log:
	       final = str(LOGS_PATH / "Mapping" / "Orfium" / "Star"/ "{sample}_Orfium_Log.final.out"),
	       log = str(LOGS_PATH / "Mapping" / "Orfium" / "Star"/ "{sample}_Orfium_Log.out"),
	       sj = str(LOGS_PATH / "Mapping" / "Orfium" / "Star"/ "{sample}_Orfium_SJ.out.tab"),
	       prog = str(LOGS_PATH / "Mapping" / "Orfium" / "Star"/ "{sample}_Orfium_Log.progess.out"),
	       star = str(LOGS_PATH / "Mapping" / "Orfium" / "Star"/ "{sample}_star.out"), 
	       bowtie2_out = str(LOGS_PATH / "Mapping" / "Orfium" / "Star"/ "{sample}_bowie2_mapping.out")
	    benchmark:
	       str(BENCHMARKS_PATH / "Mapping" / "Orfium" / "Star" / "{sample}_STAR_Bowtie2_Mapping_Orfium.benchmark.txt")
	    shell:
	       "set +eu && "
	       ". $(conda info --base)/etc/profile.d/conda.sh && "
	       "conda activate {params.conda_env} && "
	       "STAR --readFilesCommand zcat " 
	       " --outSAMstrandField intronMotif "
	       " --outReadsUnmapped Fastx "
	       " --genomeDir {input.index_star}"
	       " --runThreadN {threads} " 
	       " --readFilesIn {input.fastq} " 
	       " --outFileNamePrefix {params.prefix} 2>> {log.star}; " 
	       "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -U {output.met1} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}" 
	        
rule compressed_unmapped_ALL_RNA:
   input:
     met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1")
   output:
     met1_compressed = str(DATA_PROCESSING_PATH / "Orfium_Fastq" / "{sample}" /"{sample}_Unmapped.out.mate1.fastq.gz")
   shell:
     "gzip -c {input.met1} > {output.met1_compressed}"
