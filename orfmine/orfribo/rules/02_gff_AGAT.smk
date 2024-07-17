#### === GFF Verification === ####


rule name_CDS:
    input:
        gff = str(GFF_PATH)
    output:
        gff_namedCDS = str(DATA_PROCESSING_PATH / "Edited_Gff" / ("Named.CDS_" + Path(str(GFF_PATH)).name))
    log: 
        str(LOGS_PATH / "Edited_Gff" / "Named.CDS_agat.log") 
    shell:
        "agat_convert_sp_gxf2gxf.pl -g  {input} -o {output} &> {log} "


