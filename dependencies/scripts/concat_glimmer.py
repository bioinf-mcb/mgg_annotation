"""
Reads all of the glimmer predicted orfs and save them in one table
"""

from pathlib import Path
import pandas as pd
from io import StringIO


def load_glimmer_results(glim):
    """ Load results from glimmer for multiple (>1) input records! """

    # read file
    with open(glim, 'r+') as f:
        file = f.readlines()

    # load tables
    comments, dfs = [], []
    for i, line in enumerate(file):
        line = line.strip()

        # get comments and tables
        if line[0] == '>':
            comments.append(line)

            # save tables
            if len(comments) > 1:
                rows = StringIO('\n'.join(rows))
                df = pd.read_csv(rows, header=None, delim_whitespace=True)
                dfs.append(df)
            rows = []

        # save last table
        elif i + 1 == len(file):
            rows.append(line) # add last row
            rows = StringIO('\n'.join(rows))
            df = pd.read_csv(rows, header=None, delim_whitespace=True)
            dfs.append(df)

        # table rows
        else: rows.append(line)

    # record IDs
    phageIDs = [header.split()[0].strip('>') for header in comments]

                    ###########################
                    ####### CHECKPOINTS #######
                    ###########################

    # different number of phageIDs and tables!
    if len(phageIDs) != len(dfs):
        print('Error in concat_glimmner.py!\nNumber of phage IDs dont match number of tables! Aboring!')
        exit()

    # no genes predicted for a phage!
    for df in dfs:
        if not df.shape:
            print('Error in concat_glimmner.py!\nSome phage had no genes predicted! \nRemove that phage! Aborting!')
            exit()

    # add contig name to each table
    for phageID, df in zip(phageIDs, dfs):
        df['contigID'] = phageID

    df = pd.concat(dfs)
    return df


def orient_localization(row):
    """ Swich start and stop in the table. Make lower value start and higher stop. """

    # make start always lower number
    if int(row['start']) <= int(row['stop']):
        new_start, new_stop = row['start'], row['stop']
    else:
        new_start, new_stop = row['stop'], row['start']

    if '-' in str(row['strand']): strand = '-'
    else: strand = '+'

    return pd.Series([new_start, new_stop, strand])


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


### concat glimmer results
glimmer_results = snakemake.input[0]
glimmer_output = snakemake.output[0]

# load table
df = load_glimmer_results(glimmer_results)
df.columns = ['id', 'start', 'stop', 'strand', 'score', 'contigID']

### curate table
# sort ascendingly by start, regardless strand
df[['new_start', 'new_stop', 'new_strand']] = df.apply(orient_localization, axis=1) # reorient start and stop
df.drop(['start', 'stop', 'strand'], axis=1, inplace=True) # remove old columns
df.rename({'new_start': 'start', 'new_stop': 'stop', 'new_strand': 'strand'}, axis=1, inplace=True)
df.sort_values(['contigID', 'start'], inplace=True)

# start & stop correspond to codon start & codon stop
df[['new_start', 'new_stop']] = df.apply(biological_localization, axis=1) # reorient start and stop
df.drop(['start', 'stop'], axis=1, inplace=True) # remove columns
df.rename({'new_start': 'start', 'new_stop': 'stop'}, axis=1, inplace=True)


# save table
cols = ['start', 'stop', 'strand', 'contigID']
df[cols].to_csv(glimmer_output, sep=',', index=False)
