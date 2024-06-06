"""
Organize all results in in one table.
"""

import pandas as pd

phages_tables = snakemake.input
confident_orfs_out = snakemake.output[0]

cols = ['proteinID','start','stop','strand','contigID',
        'orf_len','protein_len','codon_stop','source',
        'protein', 'orf']

dfs = []
for table in phages_tables:
    df = pd.read_csv(table)
    dfs.append(df)

df = pd.concat(dfs)
df['protein_len'] = df['protein'].apply(len)
df['proteinID'] = df.apply(lambda row: row['contigID'] + '_' + row['proteinID'], axis=1)
df = df[cols]

df.to_csv(confident_orfs_out, index=False)
