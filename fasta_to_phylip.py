from Bio import SeqIO, AlignIO

input_handle = open("/Bmo/jyakovleva/globo/maff_cdhit_150_maxalign_60.fa", "r")
output_handle = open("/Bmo/jyakovleva/globo/maff_cdhit_150_maxalign_60.phy", "w")

alignments = AlignIO.parse(input_handle, "fasta")
AlignIO.write(alignments, output_handle, "phylip-relaxed")

output_handle.close()
input_handle.close()
