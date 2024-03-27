rule samtools_filter:
    input: 
       sam_hisat2 = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Hisat2" / "Results" / "{sample}" / "{sample}.sam"),
       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Orfium" / "Bowtie2" / "Results" / "{sample}"/ "{sample}.sam")
    output: 
       bam = str(RESULTS_PATH / "BAM" / "Orfium" / "{sample}" / "{sample}.bam")
    resources:
       mem_mb = round(MEM_MB / 3)
    threads:
       MEM_MB
    params:
       sam = str(RESULTS_PATH / "BAM" / "Orfium" / "{sample}" / "{sample}.sam")
    benchmark:
       str(BENCHMARKS_PATH / "BAM" / "Orfium" / "{sample}_bam_Orfium.benchmark.txt")
    shell:
       "set +o pipefail ;"
       "grep '^@' {input.sam_hisat2} 1> {params.sam} ;"
       " grep -v '^@' {input.sam_hisat2} | grep -v 'ZS:i:' | egrep -i 'XM:i:0|XM:i:1' 1>> {params.sam} ;"
       " grep -v '^@' {input.sam_bowtie2} | grep -v 'XS:i:' | egrep -i 'XM:i:0|XM:i:1' 1>> {params.sam} ;" 
       "samtools view -@ 20 -F 3844 -q 1 -h -b {params.sam} | samtools sort -@ 20 -o {output.bam} ;" 
       " rm {params.sam};"
       
rule samtools_index:
    input: 
       bam = str(RESULTS_PATH / "BAM" / "Orfium" / "{sample}" / "{sample}.bam") 
    output: 
       bai = str(RESULTS_PATH / "BAM" / "Orfium" / "{sample}" / "{sample}.bam.bai")
    shell: 
       "samtools index {input.bam}" 


