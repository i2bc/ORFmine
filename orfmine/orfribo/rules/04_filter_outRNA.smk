#### OUT RNA USING BOWTIE2 #####

#### === OUT RNA ===== ####


rule index_outRNA_BOWTIE2:
    input:
         fasta= str(RNA_TO_EXCLUDE_PATH)
    output:
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    log:
         str(LOGS_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index_Filter_Unwanted_Sequence.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index_Filter_Unwanted_Sequence.txt")
    params:
         str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index" / "index_bowtie2")
    threads:
         THREADS_NB
    shell:
         "bowtie2-build --threads {threads} {input.fasta} {params} &> {log} ;"




rule bowtie_run_OutRNA:
    input:
       fastq = str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")),
       index = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    output:
       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Results" / "{sample}"/ "{sample}_Unmapped.fastq.gz"),
    params:
       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index" / "index_bowtie2")
    log: 
       bowtie2_out = str(LOGS_PATH /  "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence.log"),
    benchmark:
       str(BENCHMARKS_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence.benchmark.txt")
    resources: 
       mem_mb= MEM_MB
    threads:
       THREADS_NB
    shell:
       "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -U {input.fastq} --un-gz {output.sam_bowtie2} > /dev/null 2>> {log.bowtie2_out}" 



