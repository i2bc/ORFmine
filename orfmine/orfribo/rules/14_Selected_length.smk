rule select_read_lengths:
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
                if [ -n "{params.mean}" ]; then
                    selected_length --dir {params.orfstats_dir} --mean --threshold {params.mean} --output {output.table}
                fi

                if [ -n "{params.median}" ]; then
                    selected_length --dir {params.orfstats_dir}  --median --threshold {params.median} --output {output.table}
                fi

                if [ -n "{params.median}" ] && [ -n "{params.mean}" ]; then
                    selected_length --dir {params.orfstats_dir}  --both --threshold {params.mean} --output {output.table}
                fi

        """

