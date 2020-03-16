import pandas as pd
from Bio import Entrez, SeqIO

# будущие аргументы
kaiju_path = '/home/yulia/PycharmProjects/g_paramecii/kaiju.out'
contigs_path = '/home/yulia/PycharmProjects/g_paramecii/contigs500bp.fasta'
out_path = '/home/yulia/PycharmProjects/g_paramecii/unclassified_contigs.fasta'

# парсим имена неклассифицированных контихов
kaiju_out = pd.read_csv(kaiju_path, sep='\t', header=None)
nodes = list(set(kaiju_out.loc[kaiju_out[0] == 'U'][1]))

# записываем контиги в отдельный файл
contigs = SeqIO.parse(contigs_path, "fasta")
with open(out_path, 'w') as out_f:
    for seq in contigs:
        if seq.name in set(nodes):
            SeqIO.write(seq, out_f, "fasta")
