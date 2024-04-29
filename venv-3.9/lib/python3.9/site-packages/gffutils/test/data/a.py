import gffutils

db = gffutils.create_db('issue_197.gff', ':memory:', merge_strategy='error')
genes = list(db.features_of_type('gene'))

genes = list(db.merge(genes))

igss = list( db.interfeatures(genes,new_featuretype='intergenic_space') )

def transform(f):
    f['ID'] = [ '-'.join(f.attributes['ID']) ]
    return f

print('------')
for i in igss:
    print(transform(i))
print('------')

db = db.update(igss, transform=transform, merge_strategy='error')

for i in db.all_features(order_by=('seqid', 'start')):
    print(i)
