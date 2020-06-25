import pandas as pd
desired_width=320
pd.set_option('display.width', desired_width)
pd.set_option('display.max_columns', 10)
import urllib
import os
from Bio import Entrez
Entrez.email = "yakovleva.spbu@gmail.com"

type='gff'
prefix='_genomic'
acc_list = ['NZ_LS483319.1', 'NZ_CP043302.1	', 'NZ_CP026961.1']


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


for acc in acc_list:
    for id in get_ids(acc):
        summary = get_raw_assembly_summary(id)
        url = summary['FtpPath_RefSeq']
        label = os.path.basename(url)
        link = ''.join(f'{url}/{label}{prefix}.{type}.gz')
        file_path = f'{label}.{type}.gz'
        print(f'Download {label}')
        urllib.request.urlretrieve(link, file_path)