"""
Reads all of the phanotate predicted orfs and save them in one table
"""

from pathlib import Path
import pandas as pd

def orient_localization(row):
    """ Swich start and stop in the table. Make lower value start and higher stop. """

    # make start always lower number
    if int(row['start']) <= int(row['stop']):
        new_start, new_stop = row['start'], row['stop']
    else:
        new_start, new_stop = row['stop'], row['start']

    return pd.Series([new_start, new_stop])


def biological_localization(row):
    """ Switch start and stop in the table accordinlgy to strand. Stop is always codon stop, start is always codon start. """

    # make start and stop locations correspornd to start & stop codons

    # strand +; start is lower number; stop is higher number
    if row['strand'] == '+':
        # start is always lower number
        if int(row['start'] <= row['stop']): new_start, new_stop = row['start'], row['stop']
        else: new_start, new_stop = row['stop'], row['start']

    # strand -; start is higher number; stop is lower number
    elif row['strand'] == '-':
        if int(row['start'] >= row['stop']): new_start, new_stop = row['start'], row['stop']
        else: new_start, new_stop = row['stop'], row['start']

    # Something is wrong.
    else:
        row['strand'], row['start'], row['stop'] = 'WARNING! Unrecognized character.', 0, 0

    return pd.Series([new_start, new_stop])

### concat phanotate results
phanotate_results = snakemake.input
phanotate_output = snakemake.output[0]
cols = ['start', 'stop', 'strand', 'contigID']

# load table
dfs = []
for phanot in phanotate_results:
    try:
        df = pd.read_csv(phanot, sep='\t', comment='#', header=None)
        dfs.append(df)
    except pd.errors.EmptyDataError:
        print(f'Empty dataframe in phanotate! {path}')
        print('Try to rerun this file (remove phanotate file and rerun).')
        print('If it does not work talk to Janusz : D')


# curate table
df = pd.concat(dfs) # concatenate data strands
df.columns = ['start', 'stop', 'strand', 'contigID', 'score', 'empty'] # name columns
df[['new_start', 'new_stop']] = df.apply(orient_localization, axis=1) # reorient start and stop

# sort ascendingly by start, regardless strand
df.drop(['empty', 'score', 'start', 'stop'], axis=1, inplace=True) # remove columns
df.rename({'new_start': 'start', 'new_stop': 'stop'}, axis=1, inplace=True)
df.sort_values(['contigID', 'start'], inplace=True)

# start & stop correspond to codon start & codon stop
df[['new_start', 'new_stop']] = df.apply(biological_localization, axis=1) # reorient start and stop
df.drop(['start', 'stop'], axis=1, inplace=True) # remove columns
df.rename({'new_start': 'start', 'new_stop': 'stop'}, axis=1, inplace=True)


# save table
df[cols].to_csv(phanotate_output, sep=',', index=False)
