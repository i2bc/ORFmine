rule ORFget:
    input:
        fasta = str(FASTA_PATH),
        gff = str(DATA_PROCESSING_PATH / "Edited_Gff" / ("Named.CDS_" + Path(str(GFF_PATH)).name))
    output:
        fasta = str(DATA_PROCESSING_PATH / "Exome" / f"Exome_elongated.nfasta"),
        gff = str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.gff"),
        gff_with_genes = str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated_with_gene_features.gff")
    params:
        path = str(DATA_PROCESSING_PATH / "Exome"),
        outname = "Exome_elongated",
        features = GFF_ELEMENT_TO_COUNT
    log:
        orf_get = str(LOGS_PATH / "Exome"/ "Exome_ORFget.log")
    benchmark:
        str(BENCHMARKS_PATH / "Exome" / "Exome_ORFget.benchmark.txt")
    shell:
        "gff2prot --fna {input.fasta} --gff {input.gff} --features {params.features} --outdir {params.path} --out-basename {params.outname} --nucleic --elongate 50 --stop-end; "
        "mv {output.gff} {output.gff_with_genes};"
        "awk '$3 !~ /gene/' {output.gff_with_genes} > {output.gff};"



rule Exome_construction_gtf:
    input:
        fasta = str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.nfasta"),
        gff_with_genes = str(DATA_PROCESSING_PATH / "Exome"/ "Exome_elongated_with_gene_features.gff")
    output:
        fasta = str(DATA_PROCESSING_PATH / "Exome" / ("Exome_elongated.exons_" + Path(str(FASTA_PATH)).name)),
        gtf = str(DATA_PROCESSING_PATH / "Exome" / ("Exome_elongated.exons_" +  Path(str(GFF_PATH)).stem + ".gtf"))
    params:
        tmp_gtf = str(DATA_PROCESSING_PATH / "Exome" / "Exome_elongated.exons_tmp.gtf")
    log:
        samtools_index = str(LOGS_PATH / "Exome"/ "Exome_samtools_index.log"),
        gffread_gtf = str(LOGS_PATH / "Exome" / "Exome_gffread.log"),
        sed = str(LOGS_PATH /"Exome" / "Exome_sed.log"),
        awk = str(LOGS_PATH / "Exome" / "Exome_awk.log")
    benchmark:
        str(BENCHMARKS_PATH / "Benchmark" /  "Exome" / "Exome.benchmark.txt")
    shell:
        "samtools faidx {input.fasta} 2> {log.samtools_index} ;"
        "gffread -F -T -w {output.fasta} -o {params.tmp_gtf} -g {input.fasta} {input.gff_with_genes} 2> {log.gffread_gtf} ;"
        "sed -i 's/description[^\;]*\;//' {params.tmp_gtf} 2>> {log.sed} ;"
        "sed -i 's/\\t[A-Z]*[_]*gene_segment\\t/\\ttranscript\\t/' {params.tmp_gtf} 2>> {log.sed} ;"
        """awk -F '\\t' '{{if(NF<=9) {{print($0);}} else {{for(field=1;field<9;field++) {{printf("%s\\t",$field);}} for(field=9;field<=NF;field++) {{printf("%s ",$field);}} printf("\\n");}}}}' {params.tmp_gtf} > {output.gtf} 2>> {log.awk} ;"""
        "rm -f {params.tmp_gtf};"
        "sed -i 's/\\s$//' {output.gtf} 2>> {log.sed} ;"
