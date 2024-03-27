
# Calls Bam2Reads functions
rule Bam2Reads_Exome: 
   input:
      gff= str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.gff"),
      bam = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam"),
      bai = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam.bai"),
      psite_table = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample}" / "psite_offset.csv") 
   output: 
      reads = str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}"/ ("{sample}_{length}/Exome" + FRAG_LENGTH_L + "_reads.tab"))
   log: 
      bam2read= str(LOGS_PATH / "Bam2Reads_Exome"  / "{sample}.{length}.bam2read.log"),
      offset_grep= str(LOGS_PATH / "Bam2Reads_Exome" / "{sample}.{length}.offset_grep.log")
   resources:
        mem_mb = round(MEM_MB / 2)
   params:
        outdir = str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}" / "{sample}_{length} /"),
        outname = "Exome" + FRAG_LENGTH_L + " -features_include " + GFF_ELEMENT_TO_COUNT,
        sample_name = "{sample}",
        read_length = "{length}"
   log:
        bam2read = str(LOGS_PATH / "Bam2Reads_Exome" / "{sample}.{length}.bam2read.log"),
        offset_grep = str(LOGS_PATH / "Bam2Reads_Exome" / "{sample}.{length}.offset_grep.log")
   benchmark:
        str(BENCHMARKS_PATH / "Bam2Reads_Exome" / "{sample}.{length}.bam2read.benchmark.txt")
   shell:
        "set +o pipefail ;"
        "offset=$(grep {params.sample_name} {input.psite_table} 2> {log.offset_grep} | grep ^{params.read_length} 2>> {log.offset_grep} | cut -f7 2>> {log.offset_grep}) 2>> {log.offset_grep} ; "
        "if [ $offset = '' ]; then offset=12 ; fi ; "
        "bam2reads -shift ${{offset}} -kmer {params.read_length} -gff {input.gff} -bam {input.bam} -outpath {params.outdir} -outname {params.outname}" + " 2> {log.bam2read} ; "
        
