""" save each protein cluster MSA to seperate file """

from pathlib import Path
from Bio import SeqIO
import pandas as pd

# paths
# work_dir = Path('/Users/januszkoszucki/MGG Dropbox/Janusz Koszucki/code/repos/mgg_annotation/test/output/3_ANNOTATION')
# PCs2proteins = Path(work_dir, 'PCs2proteins.tsv')
# msa = Path(work_dir, 'msa.a3m')

# msa4view_dir = 
# msa4search_dir = 


# paths
msa = snakemake.input.msa
PCs2proteins = snakemake.input.PCs2proteins

msa_view_dir = snakemake.output.msa4view
msa_search_dir = snakemake.output.msa4search

# create folders
Path(msa_view_dir).mkdir(exist_ok=True, parents=True)
Path(msa_search_dir).mkdir(exist_ok=True, parents=True)

# load
PCs2proteins_df = pd.read_csv(PCs2proteins, sep='\t')
records = list(SeqIO.parse(msa, 'fasta'))

# get MSA for specific protein cluster
for PC, group in PCs2proteins_df.groupby('PC'):

    msa_view = Path(msa_view_dir, f'{PC}.a3m')
    msa_search = Path(msa_search_dir, f'{PC}.a3m')
    proteinIDs = group['proteinID'].to_list()

    PC_records = []
    for record in records:
        if record.id in proteinIDs:
            PC_records.append(record)
        if len(PC_records) == len(proteinIDs): break


    # save
    SeqIO.write(PC_records, msa_view, 'fasta')

    # make sure that hhsuite output is nice
    first, rest = PC_records[0], PC_records[1:]

    first.id, first.name, first.description = PC, '', ''
    PC_records = [first] + rest

    # save
    SeqIO.write(PC_records, msa_search, 'fasta')

