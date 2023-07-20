"""
Reads all of the prodigal predicted orfs and save them in one table.

Will encounter errors if (look checkpoints):
- no genes were predicted for genome (all genomes have to have at least one gene)
- no header was found in comments from prodigal table
"""

from pathlib import Path
import pandas as pd
import re

def load_prodigal_results(prod):
    """ Load results from prodigal for multiple (>1) input records! """

    # read file
    with open(prod, 'r+') as f:
        file = f.readlines()

    # load tables
    comments, dfs = [], []
    for i, line in enumerate(file):
        line = line.strip()

        # get comments and tables
        if line[0] == '#':
            comments.append(line)

            # save tables
            if len(comments) > 2 and len(comments) % 2 != 0:
                df = pd.DataFrame({'loc': rows})
                dfs.append(df)
            rows = []

        # save last table
        elif i + 1 == len(file):
            rows.append(line) # add last row
            df = pd.DataFrame({'loc': rows})
            dfs.append(df)

        # table rows
        else: rows.append(line)

    # record IDs
    headers = [re.search("\".*\"", header).group() for header in comments[::2]]
    phageIDs = [header.strip('\"').split()[0] for header in headers]

                    ###########################
                    ####### CHECKPOINTS #######
                    ###########################

    # different number of phageIDs and tables!
    if len(phageIDs) != len(dfs):
        print('Error in concat_prodigal.py!\nNumber of phage IDs dont match number of tables! Aboring!')
        exit()

    # no genes predicted for a phage!
    for df in dfs:
        if not df.shape:
            print('Error in concat_prodigal.py!\nSome phage had no genes predicted! \nRemove that phage! Aborting!')
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


# paths
prodigal_results = snakemake.input[0]
prodigal_output = snakemake.output[0]

# load table
df = load_prodigal_results(prodigal_results)
df.columns = ['loc', 'contigID']

# curate table
df[['id', 'start', 'stop', 'strand']] = df['loc'].str.split('_', expand=True) # split strig to columns
df['start'] = pd.to_numeric(df['start'], downcast='integer')
df['stop'] = pd.to_numeric(df['stop'], downcast='integer')
df[['new_start', 'new_stop']] = df.apply(orient_localization, axis=1) # reorient start and stop

# sort ascendingly by start, regardless strand
df.drop(['start', 'stop'], axis=1, inplace=True) # remove columns
df.rename({'new_start': 'start', 'new_stop': 'stop'}, axis=1, inplace=True)
df.sort_values(['contigID', 'start'], inplace=True)

# start & stop correspond to codon start & codon stop
df[['new_start', 'new_stop']] = df.apply(biological_localization, axis=1) # reorient start and stop
df.drop(['start', 'stop'], axis=1, inplace=True) # remove columns
df.rename({'new_start': 'start', 'new_stop': 'stop'}, axis=1, inplace=True)

# save table
cols = ['start', 'stop', 'strand', 'contigID']
df[cols].to_csv(prodigal_output, sep=',', index=False)
