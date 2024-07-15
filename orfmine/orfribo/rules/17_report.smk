if config.get('rna_to_exclude'):
        rule report_stats:
             input: 
                  brut = expand(str(FASTQ_PATH /"{sample}.fastq.gz"), sample=SAMPLES),
                  trimming = expand(str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")), sample=SAMPLES),
                  #star = expand(str(DATA_PROCESSING_PATH / "Trimmed_Filtred_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz"), sample=SAMPLES),
                  bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}"/ "{sample}_Unmapped.fastq.gz"), sample=SAMPLES)
             output: 
                  str(RESULTS_PATH / "report_analysis.txt")
             params :
                  script = "scripts/summary.sh",
                  tool = MAPPING_TOOL,
                  #star = expand(str(DATA_PROCESSING_PATH / "Trimmed_Filtred_Fastq" / "{sample}" / "{sample}_Unmapped.out.mate1.fastq.gz"), sample=SAMPLES),
                  bowtie2 = expand(str(DATA_PROCESSING_PATH / "Mapping" / "Mapping_Unwanted_Sequence_And_Filtering" / "Results" / "{sample}"/ "{sample}_Unmapped.fastq.gz"), sample=SAMPLES)
             shell:
                  """
                  if [ "{params.tool}" == "hisat2" ]; then
                      sh {params.script} {output} {input.brut} {input.trimming} {params.bowtie2} 
                  elif [ "{params.tool}" == "star" ]; then
                      sh {params.script} {output} {input.brut} {input.trimming} {params.bowtie2}
                  else
                      echo "Unsupported tool: {params.tool}" >&2
                      exit 1
                  fi
                 """                       

else:
        rule report_stats: 
             input:
                  brut = expand(str(FASTQ_PATH /"{sample}.fastq.gz"), sample=SAMPLES),
                  trimming = expand(str(DATA_PROCESSING_PATH / "Trimming" / "Trimmed_fastq" / "{sample}" / ("{sample}.cutadapt" + FRAG_LENGTH_L + ".fastq.gz")), sample=SAMPLES),
             output:
                  str(RESULTS_PATH / "report_analysis.txt")
             params :
                  script = "scripts/summary.sh",
                  tool = MAPPING_TOOL
             shell:
                  " sh {params.script} {output} {input.brut} {input.trimming}"


