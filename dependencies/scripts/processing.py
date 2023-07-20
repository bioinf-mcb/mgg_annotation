"""
Script processing results from phanotate, prodigal and glimmer.

1. Core - ORFs that were identically predicted by three tools.
2. Subcore - ORFs detected identically by two tools.
3. Single - ORF detected by one tool and do not sharing codon stop with any other ORF.
4. Core-elongated - ORFs detected by three tools, sharing codon stop, but at least one ORF has codon start detected differently. The codon start maximizing length of ORF is taken.
5. Subcore-elongated - ORFs predicted by two tool which share codon stop, but do not share codon start. The codon start maximizing length of ORF is taken.

6. Single-mixed - ORFs that share codon stop, but are on different strands (probably errors?)
Remove ambigous ORFs - ORFs that consist more than 0.65 length of the total contig length.

"""


# load modules
import pandas as pd


def get_orf_len(row):
    """ calculate ORF length based on start and stop """

    if row['strand'] == '+': length = row['stop'] - row['start'] + 1
    elif row['strand'] == '-': length = row['start'] - row['stop'] + 1

    return pd.Series([length])

# paths
input_metadata = snakemake.input[0]
PHANOTATE_PATH = snakemake.input[1]
PRODIGAL_PATH = snakemake.input[2]
GLIMMER_PATH = snakemake.input[3]

# orfs = snakemake.output[0]
# common_core = snakemake.output[1]
# overlapping = snakemake.output[2]
ambigous = snakemake.output[0]
confident = snakemake.output[1]

### load tables
phanotate_df = pd.read_csv(PHANOTATE_PATH) # loading phanotate file
prodigal_df = pd.read_csv(PRODIGAL_PATH) # loading prodigal file
glimmer_df = pd.read_csv(GLIMMER_PATH) # loading glimmer

phanotate_df = phanotate_df.assign(source = lambda x: 'phanotate')
prodigal_df = prodigal_df.assign(source = lambda x: 'prodigal')
glimmer_df = glimmer_df.assign(source = lambda x: 'glimmer')


df = pd.concat([phanotate_df, prodigal_df, glimmer_df])

### process results
final_orfs = {}
for cid, contig in df.groupby('contigID'):
    for oid, orf in contig.groupby('stop'):
        # singlet detections
        if len(orf) == 1:
            source_tab = [ True if t in orf['source'].unique() else False for t in ['phanotate', 'prodigal', 'glimmer']]
            status = 'single'
            row = [orf['start'].iloc[0], oid, orf['strand'].iloc[0], cid] + source_tab + [status]
            final_orfs[str(cid) + '-' + str(oid)] = row

        # nonsingle detections (>1)
        else:
            if orf['source'].nunique() == 3:
                # common core case: identical orfs by 3 tools
                if (orf['strand'].nunique() == 1) & (orf['start'].nunique() == 1) & (orf['stop'].nunique() == 1) :
                    source_tab = [True, True, True]
                    status = 'core'
                    row = [orf['start'].iloc[0], oid, orf['strand'].iloc[0], cid] + source_tab + [status]
                    final_orfs[str(cid) + '-' + str(oid)] = row

                # elongated ORF case: 3 predictions with same stop, take longest start
                # "+" strand case
                elif (orf['strand'].nunique() == 1) & (orf['strand'].iloc[0] == '+'):
                    source_tab = [ True if t in orf['source'].unique() else False for t in ['phanotate', 'prodigal', 'glimmer']]
                    status = 'core-elongated'
                    row = [orf['start'].min(), orf['stop'].iloc[0], '+', cid] + source_tab + [status]
                    final_orfs[str(cid) + '-' + str(oid)] = row

                # elongated ORF case: 3 predictions with same stop, take longest start
                # "-" strand case
                elif (orf['strand'].nunique() == 1) & (orf['strand'].iloc[0] == '-'):
                    source_tab = [ True if t in orf['source'].unique() else False for t in ['phanotate', 'prodigal', 'glimmer']]
                    status = 'core-elongated'
                    row = [orf['start'].max(), orf['stop'].iloc[0], '-', cid] + source_tab + [status]
                    final_orfs[str(cid) + '-' + str(oid)] = row

                # something wrong, mixed "+" and "-" strands
                # singleton mixed
                elif (orf['strand'].nunique() > 1):
                    print('WARNING: mixed strands!', cid)
                    for idx, row in orf.iterrows():
                        source_tab = [ True if t in list(row['source']) else False for t in ['phanotate', 'prodigal', 'glimmer']]
                        status = 'single-mixed'
                        row = [orf['start'], orf['stop'], row['strand'], cid] + source_tab + [status]
                        final_orfs[str(cid) + '-' + str(oid)] = row

            elif orf['source'].nunique() == 2:
                # subcore case: identical orfs by 2 tools
                if (orf['strand'].nunique() == 1) & (orf['start'].nunique() == 1) & (orf['stop'].nunique() == 1):
                    source_tab = [ True if t in orf['source'].unique() else False for t in ['phanotate', 'prodigal', 'glimmer']]
                    status = 'subcore'
                    row = [orf['start'].iloc[0], oid, orf['strand'].iloc[0], cid] + source_tab + [status]
                    final_orfs[str(cid) + '-' + str(oid)] = row

                # subcore-elongated: 2 orfs on '+' strand, sharing stop, but different start
                elif (orf['strand'].nunique() == 1) & (orf['strand'].iloc[0] == '+'):
                    source_tab = [ True if t in orf['source'].unique() else False for t in ['phanotate', 'prodigal', 'glimmer']]
                    status = 'subcore-elongated'
                    row = [orf['start'].min(), orf['stop'].iloc[0], '+', cid] + source_tab + [status]
                    final_orfs[str(cid) + '-' + str(oid)] = row

                # subcore-elongated: 2 orfs on '-' strand, sharing stop, but different start
                elif (orf['strand'].nunique() == 1) & (orf['strand'].iloc[0] == '-'):
                    source_tab = [ True if t in orf['source'].unique() else False for t in ['phanotate', 'prodigal', 'glimmer']]
                    status = 'subcore-elongated'
                    row = [orf['start'].max(), orf['stop'].iloc[0], '-', cid] + source_tab + [status]
                    final_orfs[str(cid) + '-' + str(oid)] = row

                elif (orf['strand'].nunique() > 1):
                    print('WARNING: mixed strands!', cid)
                    for idx, row in orf.iterrows():
                        source_tab = [ True if t in list(row['source']) else False for t in ['phanotate', 'prodigal', 'glimmer']]
                        status = 'single-mixed'
                        row = [orf['start'], orf['stop'], row['strand'], cid] + source_tab + [status]
                        final_orfs[str(cid) + '-' + str(oid)] = row


confident_orfs = pd.DataFrame.from_dict(final_orfs, orient='index',
                                        columns=['start', 'stop', 'strand', 'contigID',
                                        'phanotate', 'prodigal', 'glimmer', 'status'])



################################################
############# Remove ambigous ORFs #############
################################################

input_df = pd.read_csv(input_metadata, sep='\t')
input_df = input_df[['contigID', 'contig_len [bp]']]

confident_orfs = confident_orfs.merge(input_df, on='contigID', how='left')
confident_orfs['half_contig_len [bp]'] = confident_orfs['contig_len [bp]'] * 0.65
confident_orfs['orf_len'] = confident_orfs.apply(get_orf_len, axis=1)
confident_orfs['ambigous'] = confident_orfs['orf_len'] >= confident_orfs['half_contig_len [bp]']

# save tables
confident_orfs.loc[confident_orfs['ambigous']].to_csv(ambigous, index=False)
confident_orfs.loc[~confident_orfs['ambigous']].to_csv(confident, index=False)
