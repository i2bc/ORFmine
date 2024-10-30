rule select_read_lengths:
    input:
        expand(str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}" / "{sample}_{length}" / ("Exome_{length}_reads.stats")), sample=SAMPLES, length=LENGTHS)
    output:
        table = str(DATA_PROCESSING_PATH / "Selected_length" / "{sample}" / "Selected_length.txt")
    params:
        orfstats_dir = str(DATA_PROCESSING_PATH / "Bam2Reads_Exome" / "{sample}"),
        mean = ORFSTATS_THRESHOLD_MEAN,
        median = ORFSTATS_THRESHOLD_MEDIAN
    threads:
        THREADS_NB
    shell:
        """
        set +o pipefail
        if [ -n "{params.mean}" ] && [ -n "{params.median}" ]; then
            selected_length --directory {params.orfstats_dir} --both --threshold {params.mean} --output {output.table}
        elif [ -n "{params.mean}" ]; then
            selected_length --directory {params.orfstats_dir} --mean --threshold {params.mean} --output {output.table}
        elif [ -n "{params.median}" ]; then
            selected_length --directory {params.orfstats_dir} --median --threshold {params.median} --output {output.table}
        fi

        # Check if the output file is empty
        #if [ ! -s {output.table} ]; then
        #    echo "ERROR: No read length was selected. Please check the chosen mean/median value."
        #    exit 1
        #fi
        """

