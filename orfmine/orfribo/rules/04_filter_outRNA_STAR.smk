

#### === OUT RNA ===== ####

rule index_fasta_OutRNA:
    input:
        str(RNA_TO_EXCLUDE_PATH)
    params:
        conda_env = "star"
    output:
        directory(str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index")
    log:
        str(LOGS_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index_Filter_Unwanted_Sequence.log")
    benchmark:
        str(BENCHMARKS_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index_Filter_Unwanted_Sequence.txt")
    shell:
        "set +eu && "
        ". $(conda info --base)/etc/profile.d/conda.sh && "
        "conda activate {params.conda_env} && "
        "mkdir {output} && "
        "STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} --genomeSAindexNbases 4"

rule mapping_star_OutRNA:
    input:
       fastq = str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / "{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz"),
       index = rules.index_fasta_OutRNA.output
    output:
       sam = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Results" / "{sample}"/ "{sample}_Aligned.out.sam"),
       met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Results" / "{sample}"/ "{sample}_Unmapped.out.mate1")
    params:
        prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Index" / "{sample}_"),
        conda_env = "star"
    log: 
       final = str(LOGS_PATH /  "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence_Log.final.out"),
       log = str(LOGS_PATH /  "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence_Log.out"),
       sj = str(LOGS_PATH /  "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence_SJ.out.tab"),
       prog = str(LOGS_PATH /  "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence_Log.progess.out")
    benchmark:
       str(BENCHMARKS_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "{sample}_Filter_Unwanted_Sequence.benchmark.txt")
    shell:
       "set +eu && "
       ". $(conda info --base)/etc/profile.d/conda.sh && "
       "conda activate {params.conda_env} && " 
       "STAR --readFilesCommand zcat "
       " --outSAMstrandField intronMotif "
       " --outReadsUnmapped Fastx "
       " --genomeDir {input.index}" 
       " --runThreadN 20 "
       " --readFilesIn {input.fastq} "
       " --outFileNamePrefix {params.prefix}" 


rule compressed_unmapped_outRNA: 
   input: 
     met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Filter_Unwanted_Sequence" / "Results" / "{sample}"/ "{sample}_Unmapped.out.mate1")
   output: 
     met1_compressed = str(DATA_PROCESSING_PATH / "Trimmed_Filtred_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz")
   shell:
     "gzip -c {input.met1} > {output.met1_compressed}"

