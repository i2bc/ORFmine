#### === GFF Verification === ####


rule name_CDS:
    input:
        gff = str(GFF_PATH)
    output:
        gff_namedCDS = str(DATA_PROCESSING_PATH / "Edited_Gff" / ("Named.CDS_" + Path(str(GFF_PATH)).name))
    run:
        name_in_gff = False
	gene_id_bool = True
        if name_in_gff == True:
            db = gffutils.create_db({input.gff}, ':memory:', merge_strategy='create_unique', keep_order=True) 
            with open(output.gff_namedCDS, 'w') as fout:
                for d in db.directives:
                    fout.write('##{0}\n'.format(d))
                for feature in db.all_features():
                    if feature.featuretype == config['gff_cds_feature'] or feature.featuretype == "exon":
                        parent = list(db.parents(feature, level=1))
                        if len(parent) > 0:
                            parent = parent[0]
                            if parent.attributes.get(config['gff_name_attribut']) and not feature.attributes.get(config['gff_name_attribut']):
                                feature.attributes[config['gff_name_attribut']] = [i.replace("mRNA","cds") for i in parent.attributes.get(config['gff_name_attribut'])]
                                feature.attributes[config['gff_name_attribut']][0] + "_name"
                            if parent.attributes.get('ID') and not feature.attributes.get('ID'):
                                feature.attributes["ID"] = parent.attributes["ID"]
                                feature.attributes['ID'] = feature.attributes['ID'][0] + "_CDS"
                    fout.write(str(feature) + '\n')
        if gene_id_bool:
            print("gene_id attributes are present")
            shell("ln -s {input.gff} {output.gff_namedCDS} ;")
        else:
            print("Missing at least some 'gene_ide' attribute in this gff")
            shell("sed -i -E 's/\\s/\\t/8' {output.gff_namedCDS} ;")


