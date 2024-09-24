#### === QUALITY CONTROL === ####

rule Quality_Control:
    input:
        str(FASTQ_PATH /"{sample}.fastq.gz")
    output: 
        temp(str(DATA_PROCESSING_PATH / "Quality_control" / "Before_Trimming" / "{sample}" / "{sample}_fastqc.zip")),
        temp(str(DATA_PROCESSING_PATH / "Quality_control" / "Before_Trimming" / "{sample}" / "{sample}_fastqc.html"))
    log:
        str(LOGS_PATH / "Quality_control" / "{sample}.log")
    benchmark:
        str(BENCHMARKS_PATH / "Quality_control" / "{sample}_benchmark.txt")
    params:
       outdir = str(DATA_PROCESSING_PATH / "Quality_control" / "Before_Trimming" / "{sample}")
    shell:
        """
	 fastqc {input} --outdir {params.outdir} 2> {log} 
	"""

rule Init_multiqc: 
    input: 
        expand(str(DATA_PROCESSING_PATH / "Quality_control" / "Before_Trimming" / "{sample}" / "{sample}_fastqc.zip"), sample=SAMPLES)
    output: 
        str(DATA_PROCESSING_PATH / "Quality_control" / "Before_Trimming" / "multiqc_report.html")
    log:
        str(LOGS_PATH / "Quality_control" / "multiqc_report_before_Trimming.log")
    benchmark:
        str(BENCHMARKS_PATH / "Quality_control" / "multiqc_report_Before_Trimming_benchmark.txt")
    params: 
        str(DATA_PROCESSING_PATH / "Quality_control" / "Before_Trimming/ " )
    shell: 
        " multiqc -f {input} -o {params} . > {log} 2>&1 "
