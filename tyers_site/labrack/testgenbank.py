from Bio import SeqIO
gb_file = "examples/rg2130_pJE411_reformatted.gbk"
gb_features = ""
for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
    # now do something with the record
    print "Name %s, %i features" % (gb_record.name, len(gb_record.features))
    print repr(gb_record.seq)
    for ind in xrange(len(gb_record.features)) :
        gb_features += repr(gb_record.features[ind])
    print gb_features



from Bio import SeqIO
input_handle = open("rg2130_pJE411_reformatted.gbk", "rU")
for record in SeqIO.parse(input_handle, "genbank") :
    print record
input_handle.close()

from Bio import SeqIO
gb_record = SeqIO.read('rg2130_pJE411_reformatted.gbk', 'gb')
for gb_feature in gb_record.features:
    gene=gb_feature.qualifiers['gene'][0]
    print gene