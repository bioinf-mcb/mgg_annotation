"""
Generate binary matrices per capsuke for pyseer.

KL types presence/absence.
Protein families presensce/absence.

"""

import pandas as pd


clusters_raw = snakemake.input.pcs_raw
main_raw = snakemake.input.main

clusters = snakemake.output.pcs
clusters_binary = snakemake.output.pcs_binary
main = snakemake.output.main


### Process clustering results

# load
clusters_df = pd.read_csv(clusters_raw, sep='\t')
clusters_df.columns = ['repr', 'proteinID']

# sort & rename
size_df = clusters_df.groupby('repr').size().sort_values(ascending=False).reset_index()
size_df.columns = ['repr', 'size']

# sort PCs
clusters_df = clusters_df.merge(size_df, on='repr', how='left')
clusters_df.sort_values(['size', 'repr'], ascending=[False, True], inplace=True)
clusters_df.drop('size', axis=1, inplace=True)

# rename PCs
reprs = size_df['repr'].to_list()
format = len(str(len(reprs)))
mapper = {repr: 'PC' + f'{i+1}'.zfill(format) for i, repr in enumerate(reprs)}

clusters_df['PC'] = clusters_df['repr'].map(mapper)
clusters_df = clusters_df[['PC', 'proteinID', 'repr']]

# save
clusters_df.to_csv(clusters, sep='\t', index=False)
clusters_df.drop('repr', axis=1, inplace=True)


### Add clusters to main table
main_df = pd.read_csv(main_raw, sep='\t')
main_df = main_df.merge(clusters_df, on='proteinID', how='left')

cols = ['PC', 'proteinID', 'prophageID', 'contigID', 'genomeID', 'KL', 'ST', 'start', 'stop', 'strand', 'protein', 'orf']
main_df[cols].to_csv(main, sep='\t', index=False)

### protein faimilies presence/absence matrix
main_df.rename({'genomeID': 'sample'}, axis=1, inplace=True)

main_df['value'] = 1
clusters_binary_df = main_df[['PC', 'sample', 'value']]
clusters_binary_df = clusters_binary_df.pivot_table(index='PC', columns='sample')['value']

# make values binary
clusters_binary_df.fillna(0, inplace=True)
for col in clusters_binary_df.columns:
    clusters_binary_df[col] = pd.to_numeric(clusters_binary_df[col], downcast='integer')
clusters_binary_df = clusters_binary_df.reset_index()

# save
clusters_binary_df.to_csv(clusters_binary, sep='\t', index=False)
