rule samtools_filter_Exome:
    input: 
       sam_hisat2 = str( DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Hisat2" / "Results"  / "{sample}" / "{sample}.sam"),
       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam"),
    output: 
       bam = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam")
    resources:
       mem_mb = round(MEM_MB / 3)
    threads:
       THREADS_NB
    params: 
       sam = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.sam")
    benchmark:
       str(BENCHMARKS_PATH / "BAM" / "Exome" / "{sample}_bam_orfeum.benchmark.txt")
    shell:
        "set +o pipefail ;"
	"grep '^@' {input.sam_hisat2} | uniq 1> {params.sam} ;"
	"grep -v '^@' {input.sam_hisat2} | egrep -i 'XM:i:0|XM:i:1' 1>> {params.sam} ;" 
	"grep -v '^@' {input.sam_bowtie2} 1>> {params.sam} ;"
	"samtools view -@ {threads} -F 3588 -h -b {params.sam} | samtools sort -@ {threads} -o {output.bam} ;"
	" rm -f {params.sam}" 

rule samtools_index_Exome:
    input: 
       bam = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam")
    output: 
       #bai = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / ("{sample}_" + FRAG_LENGTH_L + ".bam.bai"))
       bam = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam.bai")
    shell: 
       "samtools index {input.bam}" 

