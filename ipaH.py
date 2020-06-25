from Bio import SeqIO

# files = ['ipaH/GCF_001007915.1_ASM100791v1_genomic.gbff', 'ipaH/GCF_001578125.1_ASM157812v1_genomic.gbff',
#          'ipaH/GCF_000012005.1_ASM1200v1_genomic.gbff', 'ipaH/GCF_001027225.1_ASM102722v1_genomic.gbff',
#          'ipaH/GCF_003855115.1_ASM385511v1_genomic.gbff']

files = ['ipaH/GCF_001578125.1_ASM157812v1_genomic.gbff']

name = 'pdeudo'

for file in files:
    recs = [rec for rec in SeqIO.parse(file, 'genbank')]


    def get_cds_feature(recs, query):
        for rec in recs:
            for feature in rec.features:
                for key, val in feature.qualifiers.items():
                    if query in str(val):
                        # print('pseudo' in feature.qualifiers)
                        print('>' + feature.qualifiers['product'][0], rec.description,
                              feature.location, 'pseudo' in feature.qualifiers)
                        # print('>' + feature.qualifiers['gene'][0], rec.description)
                        # print(feature.extract(rec.seq))

    # get_cds_feature(recs, "E3 ubiquitin--protein ligase")
    get_cds_feature(recs, "IpaH")
    # get_cds_feature(recs, "invasion protein")