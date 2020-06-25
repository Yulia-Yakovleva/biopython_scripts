from Bio import SeqIO
import seaborn as sns
sns.set_style('darkgrid')
import matplotlib.pyplot as plt

genome = 'Escherichia_coli_C4/GCF_001900315.1_ASM190031v1_genomic.fna'
insertion_seqs = 'Escherichia_coli_C4/GCA_001900315_transposase.bed'
# genome = 'Shigella_flexneri_1a_228/chromosome.fna'
# insertion_seqs = 'Shigella_flexneri_1a_228/GCA_001578125_transposase.bed'



interval = 10000

genome_record = SeqIO.parse(genome, 'fasta')

def check_intervals(file):
    intervals = []
    for line in file:
        line = line.strip()
        first_int = int(line.split('\t')[1])
        second_int = int(line.split('\t')[2])
        if first_int < second_int:
            start = first_int
            end = second_int
        else:
            start = second_int
            end = first_int
        intervals.append((start, end))
    return intervals

def split_genome_on_intervals(record, interval):
    intervals = []
    for data in record:
        genome_size = len(data.seq)
        print(genome_size)
        for i in range(0, genome_size, interval):
            intervals.append((i+1, i + interval))
    return intervals

def count_intervals(just_one_genome_interval, list_of_is_intervals):
    passed_intervals = 0
    failed_intervals = 0
    left_bound = just_one_genome_interval[0]
    right_bound = just_one_genome_interval[1]
    for el in list_of_is_intervals:
        start = el[0]
        end = el[1]
        if start > right_bound or end < left_bound:
            failed_intervals +=1
            # print('failed')
        else:
            if start < right_bound and end > left_bound:
                passed_intervals +=1
                # print('ha-ha, classic')
            elif start < left_bound and end < right_bound:
                # print('okay, you win')
                passed_intervals +=1
            elif start < right_bound and end > right_bound:
                # print('hey-hey, NOT TODAY')
                failed_intervals += 1
    return passed_intervals, failed_intervals


with open(insertion_seqs, 'r') as is_seqs, open(genome, 'r') as genome:

    insertion_seqs_intervals = check_intervals(is_seqs)
    genome_intervals = split_genome_on_intervals(genome_record, interval)

    check_results = {}

    el_count = 0
    for el in genome_intervals:
        el_count +=1
        passed, failed = count_intervals(el, insertion_seqs_intervals)
        check_results[el_count] = passed

    for int, count in zip(genome_intervals, check_results):
        print(f"NZ_CP010121.1\t{int[0]}\t{int[1]}\t{check_results[count]}")
        # print(max(check_results.values()))

    # PLOT CODE
    # plot = sns.barplot(x=list(check_results.keys()),
    #                    y=list(check_results.values()))
    # plot.set(xlabel='Genomic interval with step ' + f"{interval} bp",
    #          ylabel='Insertion sequence counts')
    # xticklabels = plot.set_xticklabels([f"{round(value[0]/1000000, 1)}-{value[1]/1000000} Mb"
    #                                     if round(value[0] / 1000000, 1) > 0
    #                                     else f"{round(value[0]/1000000)}-{value[1]/1000000} Mb"
    #                                     for value in genome_intervals])
    # # plot.set_xticklabels(xticklabels, rotation=70, horizontalalignment='right')
    # plot.set_xticklabels(xticklabels, rotation=90)
    # plt.title(f'Distribution of ISs across the $Bacteria$ genome')
    #
    # plt.tight_layout()
    # plt.show(plot)
    # # plt.savefig('IS_distribution/Shigella_sonnei_FDAARGOS_524.pdf', fromat="pdf")