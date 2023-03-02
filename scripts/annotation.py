"""
Read HHblits sensitive search results tables, filter them and save in one table.
"""

print('Executing annotation.py!')

# modules
from csb.bio.io import HHOutputParser
import pandas as pd
import numpy as np
from pathlib import Path

from utils import bcolors
from annotation_utils import get_phrog, curate_columns, get_no_hit_frames
from annotation_utils import report_phrogs, report_alan, report_pfam, report_ecod
from annotation_utils import combine_function_and_confidence

# functions
def get_hits(tables, dbname='UNKNOWN'):
    """ parse mmseqs output """
    colnames = ['query','target','qstart','qend','qlength','tstart','tend', \
                'tlength','ident','bits','evalue','prob','pvalue', 'name']

    rows, corrupted, parser = [], [], HHOutputParser()
    for table in tables:
        try:
            parsed_table = parser.parse_file(table)
        except:
            corrupted.append(table)
            continue

        qname = parsed_table.query_name
        for hit in parsed_table:
            hit_values = [qname, hit.id, hit.qstart, hit.qend, hit.qlength,
                          hit.start, hit.end, hit.slength,
                          round(int(hit.identity)/100, 2), hit.score, hit.evalue,
                          hit.probability, hit.pvalue, hit.name]

            rows.append(hit_values)

    df = pd.DataFrame(rows, columns=colnames)

    # coverage & format
    df['prob'] = np.round(df['prob'], 3)
    df['qcov'] = np.round(np.absolute(df['qstart'] - 1 - df['qend']) / df['qlength'], 2)
    df['tcov'] = np.round(np.absolute(df['tstart'] - 1 - df['tend']) / df['tlength'], 2)
    
    df['db'] = dbname # name of database
    

    # order columns
    colfinal = ['query','target','prob','pvalue','ident','qcov','tcov','bits', \
                'qstart','qend','qlength','tstart','tend', 'tlength', 'evalue', \
                'db','name']

    return df[colfinal]


# paths & params
print('Loading paths... ', end='')
pcs_table = snakemake.input.PCS_TABLE
main_input = snakemake.input.ORFS_TABLE
phrogs_tables = list(snakemake.input.phrogs)
alan_tables = list(snakemake.input.alan)
pfam_tables = list(snakemake.input.pfam)
ecod_tables = list(snakemake.input.ecod)

PHROGS_ANNOT = snakemake.params.PHROGS_ANNOT
MGG_PHROGS_ANNOT = snakemake.params.MGG_PHROGS_ANNOT
ALAN_ANNOT = snakemake.params.ALAN_ANNOT

results = snakemake.output.results
annot = snakemake.output.annotation
main_with_functions = snakemake.output.orfs
print('Done!')

# load & merge
print('Loading hhr tables... ', end='')
phrogs_df = get_hits(phrogs_tables, dbname='PHROGS')
alan_df = get_hits(alan_tables, dbname='ALANDB')
pfam_df = get_hits(pfam_tables, dbname='PFAM')
ecod_df = get_hits(ecod_tables, dbname='ECOD')

### MAP MGG PHROGS functions and ALAN function ###

# replace PHROGS functions to CUSTOM MGG PHROGS functions.
print('load PHROGS & ALAN metadata tables... ', end='')
# load
phrogs_annot_df = pd.read_csv(PHROGS_ANNOT, sep='\t') # load PHROGS
mgg_phrogs_annot_df = pd.read_csv(MGG_PHROGS_ANNOT, sep=';') # load MGG PHROGS

# rename columns
mgg_phrogs_annot_df = mgg_phrogs_annot_df.rename(columns={'funct.orig': 'annot', 'funct.new': 'annot_mgg'})[['annot', 'annot_mgg']]

# map to CUSTOM MGG PHROGS functions
phrogs_annot_df.merge(mgg_phrogs_annot_df, on='annot', how='left') \
               .drop('annot', axis=1) \
               .rename(columns={'annot_mgg': 'annot'})

# format
phrogs_annot_df = phrogs_annot_df[['phrog', 'color', 'annot', 'category']] # order columns
phrogs_annot_df = phrogs_annot_df.fillna('unknown function')


# ALANDB
alan_annot_df = pd.read_csv(ALAN_ANNOT, sep=',')
alan_annot_df = alan_annot_df.rename(columns={'definition': 'annot', 'funct': 'category', 'profile': 'alan_profile'})
alan_annot_df = alan_annot_df.drop(['abbrev', 'hmm.name', 'family.name', 'family.name.base', 'where.it.occurs', 'is.custom', 'custom.hmm.exists'], axis=1)
print('Done!')

# add MGG & PHROGs metadata
print('Add PHROGS and ALAN metadata... ', end='')
phrogs_df['phrog'] = phrogs_df.apply(get_phrog, axis=1)
phrogs_df = phrogs_df.merge(phrogs_annot_df, on='phrog', how='left')

# add ALAN metadata
alan_df['alan_profile'] = alan_df['target']
alan_df = alan_df.merge(alan_annot_df, on='alan_profile', how='left')

# master table
print('combine tables... ', end='')
df = pd.concat([phrogs_df, alan_df, pfam_df, ecod_df], axis=0)

# format columns
df[['phrog', 'alan_profile']] = df[['phrog', 'alan_profile']].fillna('0')
df['phrog'] = df['phrog'].astype(int)
df['phrog/alan_profile'] = df.apply(curate_columns, axis=1)
df = df.drop(['phrog', 'alan_profile'], axis=1)

# sort table per PC
df['tmp'] = df['query'].str.strip('PC').astype(int)
df = df.sort_values('tmp', ascending=True).drop('tmp', axis=1)

### save resutlts table ###
print('save... ', end='')
df.to_csv(results, sep='\t', index=False)
print('Done!')


#### GET BEST HIT(S) FOR EACH DB ####
best_hits_dfs = []
for pcid, pc in df.groupby('query'):
    
    # default (no hit)
    phrogs_df, alan_df, pfam_df, ecod_df = get_no_hit_frames(pcid)
    
    for dbid, db in pc.groupby('db'):
        
        # report function
        if dbid == 'PHROGS': phrogs_df = report_phrogs(db, max_evalue=10**-3, nfunc2report=2, verbose=False)
        elif dbid == 'ALANDB': alan_df = report_alan(db, min_prob=0.95, nfunc2report=2, verbose=False)
        elif dbid == 'PFAM': pfam_df = report_pfam(db, verbose=False)
        elif dbid == 'ECOD': ecod_df = report_ecod(db, verbose=False)
        else: pass

    # save function
    best_hits_dfs.append(phrogs_df)
    best_hits_dfs.append(alan_df)
    best_hits_dfs.append(pfam_df)
    best_hits_dfs.append(ecod_df)


# combine results for each db
best_hits_df = pd.concat(best_hits_dfs).reset_index(drop=True)
# simplify table
best_hits_df = best_hits_df[['query','target', 'prob', 'qcov', 'tcov', 'bits', 'evalue', 'report_label', 'report_function', 'report_params', 'report_confidence']]
# add confidence to function string
best_hits_df['report_function'] = best_hits_df.apply(combine_function_and_confidence, axis=1)
# drop confidence
best_hits_df = best_hits_df.drop('report_confidence', axis=1)
# save
best_hits_df.to_csv(annot, sep='\t', index=False)


#### ADD FUNCTIONS TO MAIN (PROTEINS) TABLE

# read
main_df = pd.read_csv(main_input)
pcs_df = pd.read_csv(pcs_table, sep='\t')


# add protein clusters
main_df = main_df.merge(pcs_df, on='proteinID', how='left')

# reformat
best_hits_df['db_function'] = best_hits_df.apply(lambda row: row['report_label'] + '_' + 'function' ,axis=1)
best_hits_df['db_params'] = best_hits_df.apply(lambda row: row['report_label'] + '_' + 'params' ,axis=1)

# reshape tables
function_df = best_hits_df.pivot(columns='db_function', index='query', values='report_function').copy()
params_df = best_hits_df.pivot(columns='db_params', index='query', values='report_params').copy()

# rename index
function_df.index.name = 'PC'
params_df.index.name = 'PC'

# combine
report_df = function_df.merge(params_df, on='PC', how='left')

# map function
main_df = main_df.merge(report_df, on='PC', how ='left')
main_df = main_df.fillna('-')
main_df.to_csv(main_with_functions, sep='\t', index=False)

