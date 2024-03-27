# Calls ORFstat functions
rule ORFstats:
    input:
        psite_table = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample}" / "psite_offset.csv"),
        reads = str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}"/ ("{sample}_{length}/Exome" + FRAG_LENGTH_L + "_reads.tab"))
    output:
        stats = str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}" / "{sample}_{length}" / ("Exome" + FRAG_LENGTH_L + "_reads.stats"))
    log:
        orfstats = str(LOGS_PATH / "ORFstat" / "{sample}.{length}.orfstats.log")
    benchmark:
        str(BENCHMARKS_PATH / "ORFstat" / "{sample}.{length}.orfstats.txt")
    params:
        #sample_name = "{sample}",
        outdir = str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}" / "{sample}_{length}"),
        read_length = "{length}"
    shell:
        "set +o pipefail;"
        "orfstats -tab {input.reads} -N 10 -out {params.outdir} 2> {log.orfstats} ;"

