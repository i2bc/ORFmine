# Performs qualitative analysis with riboWaltz
rule riboWaltz_Exome:
    input:
        Exome_gtf = str(DATA_PROCESSING_PATH / "Exome" / ("Exome_elongated.exons_" +  Path(str(GFF_PATH)).stem + ".gtf")),
        #config= pkg_resources.resource_filename("orfribo", "config.yaml"),
        #config=config.get("config", ""),
        bam_folder = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}" / "{sample}.bam.bai")
    output:
        psite_table = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample}" / "psite_offset.csv")
    resources:
        mem_mb = MEM_MB
    priority: 11
    params: 
        bam_folder = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}"),
	psite_dir = str(RESULTS_PATH / "Psite /"),
	ribo = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample} /"),
        psite = str(RESULTS_PATH / "Psite" / "{sample}_psite_table.csv"),
        min_length = MIN_READ_LENGTH,
        max_length = MAX_READ_LENGTH
    shell:
        "touch {output.psite_table} ; "
        "Rscript {periodicity_riboWaltz_exome} {input.Exome_gtf} {params.bam_folder} {params.min_length} {params.max_length}  {params.ribo} ; "
        "rm -f {OUT_BASE_PATH}/Rplots.pdf ; "
        "mkdir -p {params.psite_dir} ; "
        "cp {output.psite_table} {params.psite}"
        
