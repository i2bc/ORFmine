#### === QUALITY CONTROL === ####

rule make_fastqc:
    input:
        str(FASTQ_PATH /"{sample}.fastq.gz")
    output:
        str(DATA_PROCESSING_PATH / "Quality_control" / "{sample}" / "{sample}_fastqc.zip"),
        str(DATA_PROCESSING_PATH / "Quality_control" / "{sample}" / "{sample}_fastqc.html")
    log:
        str(LOGS_PATH / "Quality_control" / "{sample}.log")
    benchmark:
        str(BENCHMARKS_PATH / "Quality_control" / "{sample}_benchmark.txt")
    params:
       outdir = str(DATA_PROCESSING_PATH / "Quality_control" / "{sample}")
    shell:
        """
	 fastqc {input} --outdir {params.outdir} 2> {log} 
	"""

