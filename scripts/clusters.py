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