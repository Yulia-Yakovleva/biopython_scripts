from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('darkgrid')
import pandas as pd

genome_1 = 'Escherichia_coli_C4/GCF_001900315.1_ASM190031v1_genomic.fna'
is_elements_1 = 'Escherichia_coli_C4/GCA_001900315_transposase.bed'
genome_2 = 'Shigella_flexneri_1a_228/chromosome.fna'
is_elements_2 = 'Shigella_flexneri_1a_228/GCA_001578125_transposase.bed'

interval = 5000


def data_to_dict(genome, is_elements, interval):
    record = SeqIO.parse(genome, 'fasta')
    intervals = []
    for data in record:
        genome_size = len(data.seq)
        print(genome_size)  # 4990476
        for i in range(0, genome_size, interval):
            intervals.append((i, i + interval))

    is_intervals = []
    with open(is_elements, 'r') as in_f:
        for line in in_f:
            line = line.strip()
            start = int(line.split('\t')[1])
            end = int(line.split('\t')[2])
            is_intervals.append((start, end))

    intrvls = []
    counts = []

    bnd = 0
    for bound in intervals:
        count = 0
        # print(bound[0], bound[1], 'this is genome')
            if interval[0] < bound[1] and interval[1] < bound[1]:
                # print('both intervals less then right bound')
                if interval[1] > bound[1]:
                    print('конец интервала больше левой границы, YES')
                    count += 1
                else:
                    pass
                    # print('конец интервала меньше левой границы, NO')
            else:
                pass
                # print('обе границы интервала должны быть меньше правой границы, NO')
        print(count)
        bnd += 1
        counts.append(count)
        intrvls.append(bnd)
    # print(intrvls)
    return intrvls, counts


intrvls_1, counts_1 = data_to_dict(genome_1, is_elements_1, interval)
intrvls_2, counts_2 = data_to_dict(genome_2, is_elements_2, interval)

dct_1 = pd.DataFrame()
dct_2 = pd.DataFrame()
dct_1['int'] = intrvls_1
dct_1['counts'] = counts_1
dct_1['Strain'] = ["E. coli" for _ in range(len(counts_1))]
dct_2['int'] = intrvls_2
dct_2['counts'] = counts_2
dct_2['Strain'] = ["Shigella" for _ in range(len(counts_2))]
all = pd.concat([dct_1, dct_2])
# print(all)
pd.DataFrame.to_csv(all, 'for_Lavr.csv')

# plot = sns.scatterplot(x='int', y='counts', hue='Strain', data=all)
# plot.set(xlabel='Genomic interval with step '+f"{interval} bp", ylabel='ISs counts')
# plt.show(sns)

# ax = sns.barplot(x="int", y="counts", data=all)
# plt.show(ax)

