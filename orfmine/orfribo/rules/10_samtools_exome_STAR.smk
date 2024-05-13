rule samtools_filter_Exome:
    input: 
       sam_star= str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Aligned.out.sam"),
       sam_bowtie2 = str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Bowtie2" / "Results" / "{sample}" / "{sample}.sam") 
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
        ## Pour STAR       
        ''' awk -F'\t' '/^@/ && !seen[$0]++ || $15 == "nM:i:1" || $15 == "nM:i:0"' {input.sam_star} > {params.sam} ;'''       
        # Pour Bowtie2        
        '''awk -F'\t' '!/^@/ && ($15 == "XM:i:1" || $15 == "XM:i:0")' {input.sam_bowtie2} 1>> {params.sam} ;'''        
        # Samtools filter        
        "samtools view -@ 20 -F 3588 -h -b {params.sam} | samtools sort -@ 20 -o {output.bam};"
        " rm {params.sam};" 

rule samtools_index_Exome:
    input: 
       bam = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam")
    output: 
       bai = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam.bai")
    shell: 
       "samtools index {input.bam}" 
