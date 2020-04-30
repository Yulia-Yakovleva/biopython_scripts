import pandas as pd
desired_width=320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)
import urllib
import os
from Bio import Entrez
Entrez.email = "yakovleva.spbu@gmail.com"

type='faa'
prefix='_protein'


def parse_shigella_ids(stats):
    # parse table ang get accession list
    shigella_stats = stats.loc[stats['Organism'].str.contains('Shigella')]
    acc_set = set(list(shigella_stats['AssemblyID']))
    return list(acc_set)


def parse_plasmid_proteins():
    pass

def get_ids(term):
    # finds the ids associated with the assembly
    ids = []
    handle = Entrez.esearch(db="assembly", term=term)
    record = Entrez.read(handle)
    ids.append(record["IdList"])
    return ids

#Fetch raw output
def get_raw_assembly_summary(id):
    handle = Entrez.esummary(db="assembly", id=id, report="full")
    record = Entrez.read(handle)
    return(record['DocumentSummarySet']['DocumentSummary'][0])


stats = pd.read_csv('total_stats.csv')
acc_list = parse_shigella_ids(stats)
for acc in acc_list:
    for id in get_ids(acc):
        summary = get_raw_assembly_summary(id)
        url = summary['FtpPath_RefSeq']
        label = os.path.basename(url)
        link = ''.join(f'{url}/{label}{prefix}.{type}.gz')
        file_path = f'phylo/{label}.{type}.gz'
        print(f'Download {label}')
        urllib.request.urlretrieve(link, file_path)
