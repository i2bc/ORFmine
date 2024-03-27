#### === TRIMMING === ####


# Find the adapter sequence if not set in config file
rule find_adapter_sequence:
    input:
        fastq = str(FASTQ_PATH / "{sample}.fastq.gz")
    output:
        adapter = str(DATA_PROCESSING_PATH / "Trimming" / "Adapters" / "{sample}" / "{sample}.txt")
    log:
        rscript = str(LOGS_PATH / "Trimming" / "Adapters" / "{sample}.rscript.log"),
        sed = str(LOGS_PATH / "Trimming" / "Adapters" / "{sample}.sed.log"),
        echo = str(LOGS_PATH / "Trimming" / "Adapters" / "{sample}.echo.log"),
        touch = str(LOGS_PATH / "Trimming" / "Adapters" / "{sample}.touch.log")
    params:
        str(DATA_PROCESSING_PATH / "Trimming" / "Adapters" / "{sample}" "/")
    shell:
        "touch {output.adapter} 2> {log.touch};"
        "if [ -z " + SEQUENCE_ADAPTER + " ]; then "
        "Rscript " + find_adapter_sequence + " {input.fastq} " + " {params} " + "2> {log.rscript} ;"
        "elif [ '" + ARE_ADAPTERS_TRIMMED + "' = 'no' ]; then echo " + SEQUENCE_ADAPTER + " 1> {output.adapter} 2> {log.echo};"
        "fi;"

rule adapt_trimming:
    input:
        fastq = str(FASTQ_PATH / "{sample}.fastq.gz"),
        adapt_seq = str(DATA_PROCESSING_PATH / "Trimming" / "Adapters" / "{sample}" / "{sample}.txt")
    output:
        cut_fastq = str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz"))
    resources:
        mem_mb = MEM_MB
    threads:
        THREADS_NB
    params: 
        min = MIN_READ_LENGTH,
        max = MAX_READ_LENGTH, 
        trimmed = ARE_ADAPTERS_TRIMMED
    log:
        trim_value = str(LOGS_PATH / "Trimming" / "Trimmed_fastq" / "{sample}_trim_value.log"),
        cutadapt = str(LOGS_PATH / "Trimming" / "Trimmed_fastq" / "{sample}_cutadapt.log"),
        cutadapt_out = str(LOGS_PATH / "Trimming" / "Trimmed_fastq" / "{sample}_adapt_trimming.log")
    benchmark: 
        str(BENCHMARKS_PATH / "Trimming" / "Trimmed_fastq" / "{sample}_adapt_trimming.txt")
    shell:
        """
        adapter_sequence=$(cat {input.adapt_seq})
        echo $adapter_sequence        
        if [ {params.trimmed} == 0 ]; then
            cutadapt -a $adapter_sequence --trimmed-only -e 0.125 -j {threads} --max-n=1 -m  {params.min}  -M {params.max} -o {output.cut_fastq} {input.fastq} 1>> {log.cutadapt_out} 2> {log.cutadapt}
        else
            ln -s {input.fastq} {output.cut_fastq}
            echo "Reads are already trimmed." > {log.cutadapt_out}
        fi
        """
