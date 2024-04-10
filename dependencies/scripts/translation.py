"""
Load phage sequence and table with processed results.
Extract ORF sequences and translates them to proteins.
Saves it in form of tables and fasta files.

This is done for each phage seperately.
"""

# import modules
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# define functions
def extract_orf_and_protein(row, seq):
    """
    For each row extract orf for phage sequence.
    Translate the orf to protein seq.
    If codon stop found in seq return row with message.
    """
    start, stop, strand = int(row['start']), int(row['stop']), row['strand']

    if strand == '+':
        orf = seq[start-1:stop]
    elif strand == '-':
        orf = seq[stop-1:start]
        orf = Seq(orf).reverse_complement()
    else:
        print('Wrong value of strand! Aborting! Error in translation.py' )

        ##############################################
        ##### CODON STOP WITHIN PROTEIN SEQUENCE #####
        ##############################################

    # catch rare cases of codon stop within protein sequence
    codons_stop = ['TAA', 'TAG', 'TGA']
    orf, stop = list(str(orf)[:-3]), str(orf)[-3:] # extract codon stop
    for index in range(0, len(orf), 3):
        codon = ''.join(orf[index:index+3])
        if codon in codons_stop:
            orf = ''.join(orf)
            prot = str(Seq(orf).translate())
            return pd.Series([orf, prot, stop, 'WARNING! Codon stop within protein seq!'])
    orf = ''.join(orf)

    prot, orf = str(Seq(orf).translate()), str(orf)

    return pd.Series([orf, prot, stop, 'correct'])


# paths
phage_fasta = snakemake.input[0]
orfs_table_complete = snakemake.input[1]

phage_output_table = snakemake.output[0]
erronous_table = snakemake.output[1]
orfs_out = snakemake.output[2]
proteins_out = snakemake.output[3]

# load phage fasta file
records = list(SeqIO.parse(phage_fasta, 'fasta'))

# checkpoint
if len(records) > 1:
    print('Inaccurate number of records in input phage fasta file. \nInspect check_input.py\nAborting!!!')
    exit()


#############################################
######### Get tables with sequences #########
#############################################

# load phage record
record = list(records)[0]
seq = record.seq

# load table with orfs
complete_df = pd.read_csv(orfs_table_complete)
contig_name = Path(phage_fasta).stem # phage file name

# could be improved by groupby instead of filtering
filt_phage = (complete_df['contigID'] == contig_name)
phage_df = complete_df.loc[filt_phage] # get subset of df for given phage

# extract orfs and translate to proteins
orf_prot_df = phage_df.apply(extract_orf_and_protein, args=([seq]), axis=1)
orf_prot_df.columns = ['orf', 'protein', 'codon_stop', 'error']


### merge metadata table with translated proteins & orfs
phage_df = pd.concat([phage_df, orf_prot_df], axis=1)

### save proteins with codon stop in the sequence
filt_nice_proteins = (phage_df['error'] == 'correct')
erronous_df = phage_df[~filt_nice_proteins]

# save erronous prots
if len(erronous_df): erronous_df.to_csv(erronous_table, index=False)
else: Path(erronous_table).touch()

### get nice proteins
phage_df = phage_df.loc[filt_nice_proteins]
phage_df.drop('error', inplace=True, axis=1)

# get protein numbers
phage_df.reset_index(drop=True, inplace=True)

# get ID
phage_df.index.name = 'ID'
phage_df.index = phage_df.index + 1
phage_df.reset_index(inplace=True)

# get proteinID
phage_df['proteinID'] = [f'PROTEIN_{i}' for i in phage_df['ID'].to_list()]

# save table per phage
phage_df.to_csv(phage_output_table, index=False)

#################################################
######### Convert tables to fasta files #########
#################################################

# convert rows to records (orfs)
orf_records = []
for contig, protID, orf_seq in zip(phage_df['contigID'], phage_df['ID'], phage_df['orf']):

    record = SeqRecord(seq=Seq(orf_seq),
                       id = f'{contig}_ORF_{protID}',
                       name = '',
                       description = '')

    orf_records.append(record)

SeqIO.write(orf_records, orfs_out, 'fasta')


# convert rows to records (proteins)
prot_records = []
for contig, protID, prot_seq in zip(phage_df['contigID'], phage_df['ID'], phage_df['protein']):
    record = SeqRecord(seq=Seq(prot_seq),
                       id = f'{contig}_PROTEIN_{protID}',
                       name = '',
                       description = '')

    prot_records.append(record)

SeqIO.write(prot_records, proteins_out, 'fasta')
