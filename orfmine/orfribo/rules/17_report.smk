if MAPPING_TOOL == "hisat2":
    rule report_analysis:
        input:
              fastq_brut = FASTQ_PATH,
              trimmed_fastq = expand(str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")), sample=SAMPLES),
              log_exome_bowtie2 = expand(str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "{sample}_bowtie2_mapping.log"), sample=SAMPLES),
              log_exome_hisat2= expand(str(LOGS_PATH / "Mapping" / "Exome" / "Hisat2" / "{sample}_hisat2_mapping.log"), sample=SAMPLES),
              log_genome_hisat2 = expand(str(LOGS_PATH / "Mapping" / "Genome" / "Hisat2" / "{sample}_hisat2_mapping.log"), sample=SAMPLES),
              log_genome_bowtie2 = expand(str(LOGS_PATH / "Mapping" / "Genome" / "Bowtie2" / "{sample}_bowtie2_mapping.log"), sample=SAMPLES),
              
        params :
              log_filtred_unwanted_seq = expand(str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "{sample}_Mapping_Unwanted_Sequence_And_Filtering.log"), sample=SAMPLES),
              exclude = RNA_TO_EXCLUDE_PATH,        
              project_name = PROJECT_NAME
        output:
            str(RESULTS_PATH / f"{params.project_name}_report_analysis.txt")
        shell:
            """
            if [ -f "{params.exclude}" ]; then
                report --fastq {input.fastq_brut} --fastq_trimmed {input.trimmed_fastq} --unwanted_sequence {params.log_filtred_unwanted_seq} --exome {input.log_exome_hisat2} {input.log_exome_bowtie2} --genome {input.log_genome_hisat2} {input.log_genome_bowtie2} --output {output}
            else
                report --fastq {input.fastq_brut} --fastq_trimmed {input.trimmed_fastq} --exome {input.log_exome_hisat2} {input.log_exome_bowtie2} --genome {input.log_genome_hisat2} {input.log_genome_bowtie2} --output {output}
            fi 
            """



if MAPPING_TOOL == "star":
    rule report_analysis:
        input:
              fastq_brut = FASTQ_PATH,
              trimmed_fastq = expand(str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")), sample=SAMPLES),
              log_exome_bowtie2= expand(str(LOGS_PATH / "Mapping" / "Exome" / "Bowtie2" / "Exome_{sample}_bowtie2_star_mapping.txt"), sample=SAMPLES),
              log_exome_star = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Exome" / "Star" / "Results" / "{sample}" / "{sample}_Log.final.out"), sample=SAMPLES),
              log_genome_star = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Genome" / "Star" / "Results" / "{sample}" / "{sample}_Log.final.out"), sample=SAMPLES),
              log_genome_bowtie2 = expand(str(LOGS_PATH / "Mapping" / "Genome" / "Star"/ "Results" / "{sample}" /"Genome_{sample}_bowtie2_star_mapping.txt"), sample=SAMPLES)
        params :
              log_filtred_unwanted_seq = expand(str(LOGS_PATH /  "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "{sample}_Mapping_Unwanted_Sequence_And_Filtering.log"), sample=SAMPLES),
              exclude = RNA_TO_EXCLUDE_PATH,
              project_name = PROJECT_NAME        
        output:
            str(RESULTS_PATH / f"{params.project_name}_report_analysis.txt")
        shell:
            """
            if [ -f "{params.exclude}" ]; then
                report --fastq {input.fastq_brut} --fastq_trimmed {input.trimmed_fastq} --unwanted_sequence {params.log_filtred_unwanted_seq} --exome {input.log_exome_star} {input.log_exome_bowtie2} --genome {input.log_genome_star} {input.log_genome_bowtie2} --output {output}
            else
                report --fastq {input.fastq_brut} --fastq_trimmed {input.trimmed_fastq} --exome {input.log_exome_star} {input.log_exome_bowtie2} --genome {input.log_genome_star} {input.log_genome_bowtie2} --output {output}
            fi 
            """
