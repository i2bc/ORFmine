

#### === Removing Unwanted Sequences, These sequences should be in fasta format  ===== ####


rule index_outRNA_BOWTIE2:
    input:
         fasta= str(RNA_TO_EXCLUDE_PATH)
    output:
         expand(str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    log:
         str(LOGS_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index_Mapping_Unwanted_Sequence_And_Filtering.log")
    benchmark:
         str(BENCHMARKS_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index_Mapping_Unwanted_Sequence_And_Filtering.txt")
    params:
         str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index" / "index_bowtie2")
    threads:
         THREADS_NB
    shell:
         "bowtie2-build --threads {threads} {input.fasta} {params} &> {log} ;"




rule OutRNA_Removing:
    input:
       fastq = str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")),
       index = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index" / "index_bowtie2.{extb}.bt2"),extb=BOWTIE2)
    output:
       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}"/ "{sample}_Unmapped.fastq.gz"),
    params:
       index_names_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index" / "index_bowtie2")
    log: 
       bowtie2_out = str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "{sample}_Mapping_Unwanted_Sequence_And_Filtering.log"),
    benchmark:
       str(BENCHMARKS_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "{sample}_Mapping_Unwanted_Sequence_And_Filtering.benchmark.txt")
    resources: 
       mem_mb= MEM_MB
    threads:
       THREADS_NB
    shell:
       "bowtie2 -x {params.index_names_bowtie2} --threads {threads} -U {input.fastq} --un-gz {output.sam_bowtie2} > /dev/null 2>> {log.bowtie2_out}" 



