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

# Function to calculate overlap percentage based on the length of the Glimmer gene
def calculate_overlap(glimmer_gene, prodigal_gene):
    start1, end1 = glimmer_gene
    start2, end2 = prodigal_gene
    overlap = max(0, min(end1, end2) - max(start1, start2))
    glimmer_length = end1 - start1
    overlap_percentage = (overlap / glimmer_length) * 100
    return overlap_percentage


# Function to filter Glimmer genes based on overlap with Prodigal genes
def filter_glimmer_genes(prodigal_df, glimmer_df):
    filtered_glimmer = []

    for _, glimmer_gene in glimmer_df.iterrows():        
        glimmer_strand = glimmer_gene['strand']
        if glimmer_strand == '+':
            glimmer_start = glimmer_gene['start']
            glimmer_stop = glimmer_gene['stop']
        elif glimmer_strand == '-':
            glimmer_start = glimmer_gene['stop']
            glimmer_stop = glimmer_gene['start']
        glimmer_contig = glimmer_gene['contigID']
        glimmer_interval = (glimmer_start, glimmer_stop)

        # Get Prodigal genes on the same contig
        prodigal_genes_on_contig = prodigal_df[prodigal_df['contigID'] == glimmer_contig]
        prodigal_genes_on_contig = prodigal_genes_on_contig.reset_index(drop=True)

        # Swap values for negative strand
        for i, row in prodigal_genes_on_contig.iterrows():
        if row['strand'] == '-':
            start_value = row['start']
            end_value = row['stop']
            prodigal_genes_on_contig.iloc[i, prodigal_genes_on_contig.columns.get_loc('start')] = end_value
            prodigal_genes_on_contig.iloc[i, prodigal_genes_on_contig.columns.get_loc('stop')] = start_value

        # Check overlap with each Prodigal gene
        overlaps = [
            calculate_overlap(glimmer_interval, (row['start'], row['stop']))
            for _, row in prodigal_genes_on_contig.iterrows()
        ]

        # Keep Glimmer gene if all overlaps are less than 5% of glimmer gene length
        if all(overlap < 5 for overlap in overlaps):
            filtered_glimmer.append(glimmer_gene)

    return pd.DataFrame(filtered_glimmer)

# Filter Glimmer genes
filtered_glimmer_df = filter_glimmer_genes(prodigal_df, glimmer_df)

# Merge prodigal with filtered non-overlap glimmer genes
df = pd.concat([prodigal_df, filtered_glimmer_df], ignore_index=True)

### process results
final_orfs = {}
for cid, contig in df.groupby('contigID'):
    for oid, orf in contig.groupby('stop'):
            row = [orf['start'].iloc[0], oid, orf['strand'].iloc[0], cid, orf['source'].iloc[0]]
            final_orfs[str(cid) + '-' + str(oid)] = row

confident_orfs = pd.DataFrame.from_dict(final_orfs, orient='index',
                                        columns=['start', 'stop', 'strand', 'contigID', 'source'])



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
