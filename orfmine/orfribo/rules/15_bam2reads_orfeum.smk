import subprocess 
rule Bam2Reads_Orfium:
    input:
        intergenic_gff = str(GFF_INTERGENIC_PATH),
        bam = str(RESULTS_PATH / "BAM" / "Orfium" / "{sample}" / "{sample}.bam"),
        bai = str(RESULTS_PATH / "BAM" / "Orfium" / "{sample}" / "{sample}.bam.bai"),
        psite_table = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample}" / "psite_offset.csv"), 
        table = str(DATA_PROCESSING_PATH / "Selected_length" / "{sample}" / "Selected_length.txt")
    output:
        reads_table= str(DATA_PROCESSING_PATH / "Bam2Reads_Orfium" / "{sample}"/ ("{sample}_{length}/Orfium_{length}_reads.tab"))
    log: 
      bam2read= str(LOGS_PATH / "Bam2Reads_Orfium"  / "{sample}.{length}.bam2read.log"),
      offset_grep= str(LOGS_PATH / "Bam2Reads_Orfium" / "{sample}.{length}.offset_grep.log")
    resources:
        mem_mb = MEM_MB
    benchmark: 
        str(BENCHMARKS_PATH / "Bam2Reads_Orfium" / "{sample}.{length}.bam2read.benchmark.txt")
    params:
        outdir = str(DATA_PROCESSING_PATH / "Bam2Reads_Orfium" / "{sample}" / "{sample}_{length}/"),
        sample_name = "{sample}",
        reads_length = "{length}",
        feature = FEATURES_TO_COUNT,
        scripts = "/data/work/I2BC/fadwa.elkhaddar/BIM/ORFMINE/ORFmine/orfmine/orfribo/scripts/BAM2Reads.py"
    shell:
        """
        lengths_list=$(cat {input.table})
        in_lengths_list=""

        for length in $lengths_list; do
            if [ "$length" == "{params.reads_length}" ]; then
                offset=$(grep '{params.sample_name}' {input.psite_table} 2> {log.offset_grep} | grep '^{params.reads_length}' 2>> {log.offset_grep} | cut -f7 2>> {log.offset_grep}) 2>> {log.offset_grep};
                python3 {params.scripts} -shift ${{offset}} -kmer {params.reads_length} -gff {input.intergenic_gff} -bam {input.bam} -outpath {params.outdir} -outname Orfium_{params.reads_length} -features_include {params.feature} 2> {log.bam2read};
                in_lengths_list="True"
                break
            fi
        done

        if [ -z "$in_lengths_list" ]; then
            touch {output.reads_table}
        fi
        """
