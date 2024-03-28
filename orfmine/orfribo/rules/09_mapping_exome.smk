rule index_Exome_BOWTIE2:
    input: 
         fasta = str(DATA_PROCESSING_PATH / "Exome" / f"Exome_elongated.nfasta")
    output: 
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" /"index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    log: 
         str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "Exome_index.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index" / "Exome_index.benchmark.txt")
    params: 
         str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Index"/ "index_bowtie2")
    threads:
         THREADS_NB
    shell: 
         "bowtie2-build --threads {threads} {input.fasta} {params} &> {log} ;"
	 
	 

rule index_Exome_HISAT2:
    input:
         fasta = str(DATA_PROCESSING_PATH / "Exome" / f"Exome_elongated.nfasta")
    output:
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" / "index_hisat2.{exth}.ht2"),exth=HISAT2)
    log:
         str(LOGS_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" /"Exome_index.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping"/ "Exome" / "Hisat2" / "Index"/ "Exome_index.benchmark.txt")
    params:
         str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" / "index_hisat2")
    threads:
         THREADS_NB
    shell:
         "hisat2-build --threads {threads} {input.fasta} {params} &> {log} ;"


if config.get('rna_to_exclude'):       
	rule Mapping_Exome_STAR_Bowtie2:
            input: 
               fastq= str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Results" / "{sample}" / "{sample}_Unmapped.fastq.gz"),
               index_hisat2 = expand(str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" / "index_hisat2.{exth}.ht2"),exth=HISAT2),
               index_bowtie2 = expand(str(DATA_PROCESSING_PATH/ "Mapping"/ "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
            output:
               sam_hisat2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Results"  / "{sample}" / "{sample}.sam"),
               sam_bowtie2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
               fastq_hisat2 = str( DATA_PROCESSING_PATH / "Exome_Fastq" / "{sample}" / "{sample}_Unmapped.fastq.gz")
            params:
               index_names_bowtie2 =str(DATA_PROCESSING_PATH/ "Mapping"/ "Exome" / "Bowtie2" / "Index" / "index_bowtie2"),
               index_names_hisat2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" / "index_hisat2"),
               multi_map = MULTIMAPPING
            threads: 
               THREADS_NB
            log:
               hisat2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Hisat2" / "{sample}_hisat2_mapping.log"),
               bowtie2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "{sample}_bowie2_mapping.log")
            benchmark:
               str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Hisat2"/ "{sample}_HISAT_Bowtie2_Mapping_Orfeum.benchmark.txt")
            shell:
               "hisat2 -x {params.index_names_hisat2} --threads {threads} -k {params.multi_mapp}  -U {input.fastq} --un-gz {output.fastq_hisat2} -S {output.sam_hisat2} 2>> {log.hisat2_out};"
               "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -k {params.multi_mapp} -U {output.fastq_hisat2} -S {output.sam_bowtie2} 2>> {log.bowtie2_out};" 

else: 
	rule Mapping_Exome_STAR_Bowtie2:
            input: 
               fastq= str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")),    
               index_hisat2 = expand(str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" / "index_hisat2.{exth}.ht2"),exth=HISAT2),
               index_bowtie2 = expand(str(DATA_PROCESSING_PATH/ "Mapping"/ "Exome" / "Bowtie2" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
            output:
               sam_hisat2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Results"  / "{sample}" / "{sample}.sam"),
               sam_bowtie2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
               fastq_hisat2 = str( DATA_PROCESSING_PATH / "Exome_Fastq" / "{sample}" / "{sample}_Unmapped.fastq.gz")
            params:
               index_names_bowtie2 =str(DATA_PROCESSING_PATH/ "Mapping"/ "Exome" / "Bowtie2" / "Index" / "index_bowtie2"),
               index_names_hisat2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Index" / "index_hisat2"),
               multi_map = MULTIMAPPING
            threads:
               THREADS_NB
            log:
               hisat2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Hisat2" / "{sample}_hisat2_mapping.log"),
               bowtie2_out = str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "{sample}_bowie2_mapping.log")
            benchmark:
               str(BENCHMARKS_PATH / "Mapping" / "Exome" / "Hisat2"/ "{sample}_HISAT_Bowtie2_Mapping_Orfeum.benchmark.txt")
            shell:
               "hisat2 -x {params.index_names_hisat2} --threads {threads}   -k {params.multi_mapp} -U {input.fastq} --un-gz {output.fastq_hisat2} -S {output.sam_hisat2} 2>> {log.hisat2_out};"
               "bowtie2 -x {params.index_names_bowtie2} --threads {threads}  -k {params.multi_mapp} -U {output.fastq_hisat2} -S {output.sam_bowtie2} 2>> {log.bowtie2_out};"

