rule index_Exome_BOWTIE2:
    input: 
         fasta = str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.nfasta")
    output: 
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2",extb=BOWTIE2))
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
    params:
        conda_env = "star"
    output:
        directory(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Index"))
    log:
        str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Index" / "Exome_index.log")
    benchmark:
        str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Star" / "Index" / "Exome_index.txt")
    shell:
        "set +eu && "
        ". $(conda info --base)/etc/profile.d/conda.sh && "
        "conda activate {params.conda_env} && "
        "mkdir {output} && "
        "STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} --genomeSAindexNbases 7"
 
if config.get('fasta_outRNA'):       
	rule Mapping_Exome_STAR_Bowtie2:
	    input: 
	       fastq= str(DATA_PROCESSING_PATH / "Trimmed_Filtred_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz"),    
	       index_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Index"),
	       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2",extb=BOWTIE2))
	    output:
	       sam_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
	       met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1"),
	       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam") 
	    params:
	       prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_"),
	       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2",
	       conda_env = "star" 
	    log:
	       final = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.final.out"),
	       log = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.out"),
	       sj = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_SJ.out.tab"),
	       prog = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.progess.out"),
	       star = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Exome_{sample}_star.out"), 
	       bowtie2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Exome_{sample}_bowie2_mapping.out")
	    benchmark:
	       str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_STAR_Bowtie2_Mapping_Exome.benchmark.txt")
	    shell:
	       "set +eu && "
	       ". $(conda info --base)/etc/profile.d/conda.sh && "
	       "conda activate {params.conda_env} && "
	       "STAR --readFilesCommand zcat " 
	       " --outSAMstrandField intronMotif "
	       " --outReadsUnmapped Fastx "
	       " --genomeDir {input.index_star}"
	       " --runThreadN 20 " 
	       " --readFilesIn {input.fastq} " 
	       " --outFileNamePrefix {params.prefix}" 
	       " --outFilterMultimapNmax 10 ; "
	       "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -k 10 -U {output.met1} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}"

else: 
	rule Mapping_Exome_STAR_Bowtie2:
	    input: 
	       fastq= str(FASTQ_PATH),    
	       index_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Index"),
	       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2",extb=BOWTIE2))
	    output:
               sam_star = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
               met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1"),
               sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam")
            params:
               prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" /"{sample}_"),
               index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "index_bowtie2",
               conda_env = "star" 
            log:
               final = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.final.out"),
               log = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.out"),
               sj = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_SJ.out.tab"),
               prog = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_Exome_Log.progess.out"),
               star = str(LOGS_PATH / "Mapping" / "Exome" / "Star" / "Exome_{sample}_star.out"),
               bowtie2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Exome_{sample}_bowie2_mapping.out")
            benchmark:
               str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Star" / "{sample}_STAR_Bowtie2_Mapping_Exome.benchmark.txt")
            shell:
               "set +eu && "
               ". $(conda info --base)/etc/profile.d/conda.sh && "
               "conda activate {params.conda_env} && "
               "STAR --readFilesCommand zcat "
               " --outSAMstrandField intronMotif "
               " --outReadsUnmapped Fastx "
               " --genomeDir {input.index_star}"
               " --runThreadN 20 "
               " --readFilesIn {input.fastq} "
               " --outFileNamePrefix {params.prefix}"
               " --outFilterMultimapNmax 10 ; "
               "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -k 10 -U {output.met1} -S {output.sam_bowtie2} 2>> {log.bowtie2_out}"

rule compressed_unmapped_Exome:
    input: 
      met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Unmapped.out.mate1")
    output: 
      met1_compressed = str(DATA_PROCESSING_PATH / "Exome_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz")
    shell:
      "gzip -c {input.met1} > {output.met1_compressed}"
      
