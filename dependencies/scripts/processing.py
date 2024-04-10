"""
Script processing results from prodigal-gv and glimmer.

Remove ambigous ORFs - ORFs that consist more than 0.65 length of the total contig length.

Keep only prodigal-gv ORFs for all the next steps in the pipeline. Need an update to find glimmer-exclusive ORFs and tag them.

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
PRODIGAL_PATH = snakemake.input[1]
GLIMMER_PATH = snakemake.input[2]

ambigous = snakemake.output[0]
confident = snakemake.output[1]

### load tables
prodigal_df = pd.read_csv(PRODIGAL_PATH) # loading prodigal file
glimmer_df = pd.read_csv(GLIMMER_PATH) # loading glimmer

prodigal_df = prodigal_df.assign(source = lambda x: 'prodigal')
glimmer_df = glimmer_df.assign(source = lambda x: 'glimmer')


#INCLUDE GLIMMER NON OVERLAPPING ORFS HERE and tag them in the glimmer column
#df = pd.concat([prodigal_df, glimmer_df])

#For now we only take prodigal-gv results and tag them as False in glimmer
df = prodigal_df

### process results
final_orfs = {}
for cid, contig in df.groupby('contigID'):
    for oid, orf in contig.groupby('stop'):
            glimmer = False
            row = [orf['start'].iloc[0], oid, orf['strand'].iloc[0], cid] + [glimmer]
            final_orfs[str(cid) + '-' + str(oid)] = row

confident_orfs = pd.DataFrame.from_dict(final_orfs, orient='index',
                                        columns=['start', 'stop', 'strand', 'contigID', 'glimmer'])



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
