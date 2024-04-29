# Performs qualitative analysis with riboWaltz
rule riboWaltz_Exome:
    input:
        Exome_gtf = str(DATA_PROCESSING_PATH / "Exome" / ("Exome_elongated.exons_" +  Path(str(GFF_PATH)).stem + ".gtf")),
        config= pkg_resources.resource_filename("orfribo", "config.yaml")
    output:
        psite_table = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample}" / "psite_offset.csv")
    resources:
        mem_mb = MEM_MB
    params: 
        bam_folder = str(RESULTS_PATH / "BAM" / "Exome" / "{sample}"),
	psite_dir = str(RESULTS_PATH / "Psite /"),
	ribo = str(DATA_PROCESSING_PATH / "RiboWaltz" / "{sample} /"),
        psite = str(RESULTS_PATH / "Psite" / "{sample}_psite_table.csv")
    shell:
        "touch {output.psite_table} ; "
        "Rscript {periodicity_riboWaltz_exome} {input.config} {input.Exome_gtf} {params.bam_folder} {params.ribo} ; "
        "rm -f {OUT_BASE_PATH}/Rplots.pdf ; "
        "mkdir -p {params.psite_dir} ; "
        "cp {output.psite_table} {params.psite_dir}"
        
