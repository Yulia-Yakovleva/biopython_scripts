import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "yakovelva.spbu@gmail.com"

# будущие аргументы
kaiju_path = '/home/yulia/PycharmProjects/g_paramecii/kaiju.out'
contigs_path = '/home/yulia/PycharmProjects/g_paramecii/contigs500bp.fasta'
out_path = '/home/yulia/PycharmProjects/g_paramecii/fungi_contigs.fasta'
desired_taxon = 'Fungi'

# парсим айдишники классифицированных контихов
kaiju_out = pd.read_csv(kaiju_path, sep='\t', header=None)
taxids = list(set(kaiju_out.loc[kaiju_out[0] == 'C'][2]))

# узнаем, что это за айдишники
handle = Entrez.efetch(db="Taxonomy", id=" ".join(str(x) for x in taxids), retmode="xml")
records = Entrez.read(handle)

# получаем список айдишников, которые принаждежат грибам
Fungi_ids = []
for record in records:
    if desired_taxon in record["Lineage"]:
        Fungi_ids.append(str(record["TaxId"]))

# получаем список контигов, которые принадлежат грибам
contig_names = []
for col, row in kaiju_out.iterrows():
    if str(row[2]) in set(Fungi_ids):
        contig_names.append(row[1])

# записываем грибные контиги в отдельный файл
contigs = SeqIO.parse(contigs_path, "fasta")
with open(out_path, 'w') as out_f:
    for seq in contigs:
        SeqIO.write(seq, out_f, "fasta")
