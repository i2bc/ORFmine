rule index_reference_BOWTIE2:
    input: 
         fasta= str(FASTA_PATH) 
    output: 
         expand(str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
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





rule index_reference_HISAT2:
    input:
         fasta= str(FASTA_PATH)
    output:
         expand(str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Hisat2" / "Index" / "index_hisat2.{exth}.ht2"),exth=HISAT2)
    log:
         str(LOGS_PATH / "Mapping" / "Orfium" / "Hisat2" / "Index" / "Orfium_index.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Orfium" / "Hisat2" / "Index" / "Orfium_index.benchmark.txt")
    params:
         str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Hisat2" / "Index" / "index_hisat2")
    threads:
         THREADS_NB
    shell:
         "hisat2-build --threads {threads} {input.fasta} {params} &> {log} ;"


 
if config.get('rna_to_exclude'):       
	rule Mapping_Genome_STAR_Bowtie2:
	    input: 
	       fastq = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Results" / "{sample}" / "{sample}_Unmapped.fastq.gz"),    
               index_hisat2 = expand(str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Hisat2" / "Index" / "index_hisat2.{exth}.ht2"),exth=HISAT2),
	       index_bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
	    output:
	       sam_hisat2 = str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Hisat2" / "Results" / "{sample}" / "{sample}.sam"),
	       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
	       fastq_hisat2 = str(DATA_PROCESSING_PATH / "Orfium_Fastq" / "{sample}" / "{sample}_Unmapped.fastq.gz")
	    params:
	       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Bowtie2"/ "Index" / "index_bowtie2"),
	       index_names_hisat2 = str(DATA_PROCESSING_PATH / "Mapping"/ "Orfium"/ "Hisat2" / "Index" / "index_hisat2")
	    threads: 
	       THREADS_NB
	    log:
	       hisat2_out = str(LOGS_PATH / "Mapping" / "Orfium" / "Hisat2" / "{sample}_hisat2_mapping.log"),
	       bowtie2_out = str(LOGS_PATH / "Mapping" / "Orfium" / "Bowtie2" / "{sample}_bowie2_mapping.log")
	    benchmark:
	       str(BENCHMARKS_PATH / "Mapping" / "Orfium" / "Hisat2" / "{sample}_HISAT_Bowtie2_Mapping_Orfium.benchmark.txt")
	    shell:
	       "hisat2 -x {params.index_names_hisat2} --threads {threads} -U {input.fastq} --un-gz {output.fastq_hisat2} -S {output.sam_hisat2} 2>> {log.hisat2_out};"
	       "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -U {output.fastq_hisat2} -S {output.sam_bowtie2} 2>> {log.bowtie2_out};" 
	        
