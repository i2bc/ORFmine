# Concatenate counts tables selected into one by sample
rule concatenate_count_tables_genome:
    input:
        tables_list = expand(str(DATA_PROCESSING_PATH / "Bam2Reads_Genome" / "{sample}"/ ("{sample}_{length}/Genome_{length}_reads.tab")), sample=SAMPLES, length=LENGTHS)
    output:
        table = str(RESULTS_PATH / "Genome" / "{sample}" / f"Genome{FRAG_LENGTH_L}.mean{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}_reads_concatenated.tab")
    log:
        str(LOGS_PATH / "Bam2Reads_Genome" /" concatenate_count_tables_genome " / ("{sample}.concatenate_mean" + f"{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}.log"))
    benchmark:
        str(BENCHMARKS_PATH / "Bam2Reads_Genome" / "concatenate_count_tables_genome" / ("{sample}.concatenate.benchmark_mean" + f"{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}.txt"))
    params:
        outdir = str(RESULTS_PATH / "Genome" / "{sample}"),
        outname = f"Genome{FRAG_LENGTH_L}.mean{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}",
        sample_name = "{sample}"
    run:
        selected_tables = []
        for table in input.tables_list:
            if os.stat(table).st_size > 0:
                if params.sample_name in str(table):
                    selected_tables.append(table)
        table_string = ' '.join(selected_tables)
        print(table_string)
        shell("concatenate -tables " + table_string + " -outpath {params.outdir} -outname {params.outname}")


# Concatenate all samples concatenated tables together :
rule concatenate_all_tables:
    input:
        tables_list = expand(rules.concatenate_count_tables_genome.output, sample=SAMPLES, length=LENGTHS)
    output:
        table = str(RESULTS_PATH / "Genome" / f"all_samples_Genome{FRAG_LENGTH_L}.mean{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}_reads_concatenated.tab")
    log:
        str(LOGS_PATH / "Genome" / "concatenate_all_tables" / f"concatenate_mean{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}.log")
    benchmark:
        str(BENCHMARKS_PATH / "Genome" / "concatenate_all_tables" / f"concatenate.benchmark_mean{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}.txt")
    params:
        outdir = str(RESULTS_PATH / "Genome" "/"),
        outname = f"all_samples_Genome{FRAG_LENGTH_L}.mean{ORFSTATS_THRESHOLD_MEAN}_median{ORFSTATS_THRESHOLD_MEDIAN}",
    run:
        selected_tables = []
        for table in input.tables_list:
            if os.stat(table).st_size > 0:
                selected_tables.append(table)
        table_string = ' '.join(selected_tables)
        print(table_string)
        shell("concatenate -tables " + table_string + " -outpath {params.outdir} -outname {params.outname}")

