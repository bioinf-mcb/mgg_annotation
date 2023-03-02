""" Process raw mmseqs results into table"""

import pandas as pd


clusters_raw = snakemake.input[0]
clusters = snakemake.output[0]


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

# rename PCs
reprs = size_df['repr'].to_list()
format = len(str(len(reprs)))
mapper = {repr: 'PC' + f'{i+1}'.zfill(format) for i, repr in enumerate(reprs)}

clusters_df['PC'] = clusters_df['repr'].map(mapper)
clusters_df = clusters_df[['PC', 'proteinID', 'repr']]

# make nice
clusters_df['PC_integer'] = clusters_df.apply(lambda row: int(row['PC'].strip('PC')), axis=1)
clusters_df = clusters_df.sort_values('PC_integer', ascending=True)
clusters_df = clusters_df.drop('PC_integer', axis=1)

# save
clusters_df.to_csv(clusters, sep='\t', index=False)