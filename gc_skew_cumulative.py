from Bio import SeqIO
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

genomes = ['Escherichia_coli_C4/GCF_001900315.1_ASM190031v1_genomic.fna',
           'Escherichia_coli_CFSAN029787/chromosome.fna',
           'Shigella_flexneri_1a_228/chromosome.fna',
           'Shigella_sonnei_FDAARGOS_524/chromosome.fna',
           'Shigella_boydii_ATCC_9210/GCF_001027225.1_ASM102722v1_genomic.fna',
           'Shigella_dysenteriae_Sd197/chromosome.fna']

# genomes = ['Escherichia_coli_CFSAN029787/chromosome.fna']
# genomes = ['Shigella_boydii_ATCC_9210/GCF_001027225.1_ASM102722v1_genomic.fna']
genomes = ['Shigella_flexneri_1a_228/chromosome.fna']
# genomes = ['Escherichia_coli_C4/GCF_001900315.1_ASM190031v1_genomic.fna']

def Skew(symbol):
    if symbol == 'G':
        return 1
    elif symbol == 'C':
        return -1
    else:
        return 0


def SkewArray(genome):
    n = len(genome)
    array = np.zeros(n)
    for i in range(n):
        array[i] = Skew(genome[i])
    array = np.cumsum(array)
    return array


def MinimumSkew(genome):
    n = len(genome)
    array = SkewArray(genome)
    # array = array[500:]
    answ1 = np.argwhere(array == array.min())
    array = SkewArray(genome)
    # array = array[450000:]
    answ2 = np.argwhere(array == array.max())
    return answ1.tolist(), answ2.tolist()

for genome in genomes:
    genome_record = SeqIO.parse(genome, 'fasta')

    for rec in genome_record:
        print(genome, MinimumSkew(rec))
        data = SkewArray(rec)
        print('start plotting...')
        plt.plot(data)
        plt.show()