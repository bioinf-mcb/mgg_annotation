""" save each protein cluster MSA to seperate file """

from pathlib import Path
from Bio import SeqIO
import pandas as pd

# paths
msa_all = snakemake.input[0]
msa_view = snakemake.output.msa4view
msa_search = snakemake.output.msa4search

pc = Path(msa_view).stem
pc_table = snakemake.params[0]

# load
pcs_df = pd.read_csv(pc_table, sep='\t')
records = SeqIO.parse(msa_all, 'fasta')

# get MSA for specific protein cluster
filt = (pcs_df['PC'] == pc)
proteinIDs = pcs_df.loc[filt, 'proteinID'].to_list()

PC_records = []
for record in records:
    if record.id in proteinIDs:
        PC_records.append(record)

# save
n1 = SeqIO.write(PC_records, msa_view, 'fasta')

# make sure that hhsuite output is nice
first, rest = PC_records[0], PC_records[1:]

first.id, first.name, first.description = '', '', ''
first.id = pc

PC_records = [first] + rest

# save
n2 = SeqIO.write(PC_records, msa_search, 'fasta')

# checkpoint
if len(proteinIDs) != n1 or len(proteinIDs) != n2:
    print('Error in split_msa.py')
