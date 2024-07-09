import pandas as pd
from pathlib import Path
from Bio import SeqIO
from tabulate import tabulate
import yaml


# Define colors for printing.
class bc:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def run_checkpoint(config):

    # paths & params
    PHAGES_DIR = config['PHAGES_DIR']
    OUTPUT_DIR = config['OUTPUT_DIR']
    EXTENSION = config['INPUT_EXTENSION']
    PHAGE_MIN_LENGTH = config['PHAGE_MIN_LENGTH']

    IN_DIR_PROCESSED = Path(OUTPUT_DIR, '1_PROCESSED_INPUT')

    PHROGS_DIR = Path(config['HHSUITE']['PHROGS']).parent
    PFAM_DIR = Path(config['HHSUITE']['PFAM']).parent
    ECOD_DIR = Path(config['HHSUITE']['ECOD']).parent
    ALANDB_DIR = Path(config['HHSUITE']['ALANDB']).parent

    DB_DIRS = [PHROGS_DIR, PFAM_DIR, ECOD_DIR, ALANDB_DIR] 

    ### checkoints    
    
    # databases' files
    [check_db(DB_DIR) for DB_DIR in DB_DIRS]

    # mode [fasta/genbank]
    MODE = get_mode(EXTENSION)

    # clean & copy input records [PHAGES_DIR -> clean -> IN_DIR_PROCESSED]
    process_records(PHAGES_DIR, IN_DIR_PROCESSED, EXTENSION, PHAGE_MIN_LENGTH)

    return MODE



def process_records(PHAGES_DIR, IN_DIR_PROCESSED, EXTENSION, PHAGE_MIN_LENGTH):

    # file paths
    files = list(Path(PHAGES_DIR).glob(f'*{EXTENSION}'))

    # mode/file type [fasta/genbank]
    ftype = get_mode(EXTENSION)

    # prepare output directory
    Path(IN_DIR_PROCESSED).mkdir(exist_ok=True, parents=True)


    ### files/records metadata
    phages_metrics = []
    for file in files:
        metrics = record_metrics(file, ftype)
        phages_metrics.append(metrics)

    metrics_df = pd.DataFrame(phages_metrics, columns=['fname', 'recordID', 'ncontigs', 'length'])

    ### clean

    # filter
    one_contig = (metrics_df['ncontigs'] == 1)
    good_length = (metrics_df['length'] >= PHAGE_MIN_LENGTH)

    #non-unique recordID

    frq_df = metrics_df.groupby('recordID').size()
        
    unique_phages = frq_df[frq_df > 1].index
    
    # Wybierz phagi, które mają więcej niż jedno wystąpienie identyfikatora phag w metrics_df
    repeated_phages_df = metrics_df[metrics_df['recordID'].isin(unique_phages)]

    # Dołącz powtarzające się phagi do discard_phages_df


    discard_phages = ~(one_contig & good_length ) 

    # discard phages
    discard_phages_df = metrics_df.loc[discard_phages].copy()
    discard_phages_df = pd.concat([discard_phages_df, repeated_phages_df], ignore_index=True)

    # print to console
    print(f'{bc.WARNING}DiscardedPhagesWarning:. {bc.ENDC}')
    print(tabulate(discard_phages_df, headers='keys', tablefmt='grid'))

    ### phages to analyze
    phages_df = metrics_df.loc[~discard_phages].copy()


    # files
    print(f'{bc.OKGREEN}CopyCleanPhages... {bc.ENDC}' , end='')
    for row in phages_df.itertuples():

        # paths
        infile = Path(PHAGES_DIR, row.fname + f'.{EXTENSION}')
        outfile = Path(IN_DIR_PROCESSED, row.recordID + '.fasta')

        # load
        record = list(SeqIO.parse(file, ftype))[0]

        # overwrite orignal record ID by clean recordID
        record.id = row.recordID

        # save file with recordID
        SeqIO.write(record, outfile, 'fasta')
    
    print(f'{bc.OKGREEN}Done!{bc.ENDC}')



def record_metrics(file, ftype):
    """ metrics of phage records """

    # load
    records = list(SeqIO.parse(file, ftype))
    fname = Path(file).stem
    
    recordID, ncontigs, length = None, 0, 0
        
    ### metrics
    if records:
        record = records[0]
        recordID = clean_record_id(record.id)
        ncontigs = len(records)
        length = len(record.seq)

    # return
    metrics = [fname, recordID, ncontigs, length]
    
    return metrics


def clean_record_id(recordID):
    """ access record id and remove whitecharacters """
    return ''.join(recordID.split()).replace('.', '_')


def check_db(DATABASE_DIR):

    # not implemented: should verify specific files
    db_files = ['_a3m.ffdata','_a3m.ffindex','_cs219.ffdata','_cs219.ffindex','_hhm.ffdata','_hhm.ffindex']

    # checkpoint
    if not list(Path(DATABASE_DIR).glob('*')):
        print(f'{bc.WARNING}DatabaseEmptyWarning:. {DATABASE_DIR} {bc.ENDC}', end='\n')

def get_mode(EXTENSION):

    # params
    if EXTENSION == 'fasta': ftype = 'fasta'
    elif EXTENSION == 'gb': ftype = 'genbank'

    return ftype
    