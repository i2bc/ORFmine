

#### === OUT RNA ===== ####

rule index_fasta_OutRNA:
    input:
        str(RNA_TO_EXCLUDE_PATH)
    output:
        directory(str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index"))
    log:
        str(LOGS_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index_Mapping_Unwanted_Sequence_And_Filtering.log")
    benchmark:
        str(BENCHMARKS_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Index_Mapping_Unwanted_Sequence_And_Filtering.txt")
    shell:
        "mkdir {output} && "
        "../STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 20 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input} --genomeSAindexNbases 4"

rule mapping_star_OutRNA:
    input:
       fastq = str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")),
       index = rules.index_fasta_OutRNA.output
    output:
       sam = str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}"/ "{sample}_Aligned.out.sam"),
       met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}"/ "{sample}_Unmapped.out.mate1")
    params:
       prefix = str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" /"{sample}"/ "{sample}_"),
       introns = INTRONS_LENGTH 
    log: 
       final = str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" /"Results"/ "{sample}"/ "{sample}_Mapping_Unwanted_Sequence_And_Filtering_Log.final.out"),
       log = str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" /"{sample}"/ "{sample}_Mapping_Unwanted_Sequence_And_Filtering_Log.out"),
       sj = str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" /"{sample}"/ "{sample}_Mapping_Unwanted_Sequence_And_Filtering_SJ.out.tab"),
       prog = str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" /"{sample}"/ "{sample}_Mapping_Unwanted_Sequence_And_Filtering_Log.progess.out")
    benchmark:
       str(BENCHMARKS_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "{sample}_Mapping_Unwanted_Sequence_And_Filtering.benchmark.txt")
    shell:
       "../STAR-2.7.11b/bin/Linux_x86_64/STAR --readFilesCommand zcat "
       " --alignIntronMax {params.introns} "
       " --outReadsUnmapped Fastx "
       " --genomeDir {input.index}" 
       " --runThreadN 20 "
       " --readFilesIn {input.fastq} "
       " --outFileNamePrefix {params.prefix}" 


rule compressed_unmapped_outRNA: 
   input: 
     met1 = str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}"/ "{sample}_Unmapped.out.mate1")
   output: 
     met1_compressed = str(DATA_PROCESSING_PATH / "Trimmed_Filtred_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz")
   shell:
     "gzip -c {input.met1} > {output.met1_compressed}"

