import pandas as pd
from pathlib import Path
from Bio import SeqIO
import re


# Define colors for printing.
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def get_record_id(record):
    """ access record id and remove whitecharacters """
    return ''.join(record.id.split())


def checkpoint(IN_DIR, IN_DIR_PROCESSED, PHAGE_MIN_LENGTH, EXTENSION='fasta'):
    """
    Copy input files, split any multiple records to seprate files.
    Names files accordinlgy to record identifiers.
    """

    print(f'{bcolors.WARNING}Robust checkpoints need to be added... {bcolors.ENDC}', end='')

    # checkpoint
    if Path(IN_DIR_PROCESSED).exists():
        print(f'{bcolors.WARNING}Directory {IN_DIR_PROCESSED} exists! I am not processing the input phages! {bcolors.ENDC}', end='')
        return ''

    # output dir
    Path(IN_DIR_PROCESSED).mkdir(exist_ok=True, parents=True)

    no_records = []
    for phage_path in Path(IN_DIR).glob(f'*.{EXTENSION}'):

        # get records
        records = list(SeqIO.parse(phage_path, 'fasta'))
        # more then one record in file
        if len(records) > 1:
            for i, record in enumerate(records):

                # filter records by length
                if len(record.seq) <= int(PHAGE_MIN_LENGTH):
                    no_records.append(f'{phage_path} ({recordID})')
                    continue
                else: pass

                recordID = get_record_id(record) # get record identifier
                output = Path(IN_DIR_PROCESSED, f'{recordID}.{EXTENSION}') # output path
                try: record.description = record.description.replace(recordID, '') # fix description
                except: pass

                record.id = f'{recordID}' # add cleaned identifier
                n = SeqIO.write(record, output, 'fasta') # save

        # one record in file
        elif len(records) == 1:

            record = records[0] # get record
            recordID = record.id

            # filter records by length
            if len(record.seq) <= int(PHAGE_MIN_LENGTH):
                no_records.append(f'{phage_path} ({recordID})')
                continue
            else: pass


            recordID = get_record_id(record) # get record identifier
            output = Path(IN_DIR_PROCESSED, f'{recordID}.{EXTENSION}') # output path

            try: record.description = record.description.replace(recordID, '') # fix descriptions
            except: pass

            record.id = f'{recordID}' # add cleaned identifier
            n = SeqIO.write(record, output, 'fasta') # save

        # record not found
        else:
            no_records.append(f'{phage_path} ({recordID})')

    discarded = ' '.join(no_records)
    print(f'{bcolors.WARNING}Discarded phages: {discarded}{bcolors.ENDC}', end='\n')
    print(f'{bcolors.OKGREEN}Done!{bcolors.ENDC}')

    return no_records


def get_n_iterations(SEARCH_TOOL):

    if SEARCH_TOOL == 'hhblits':
        PHROGs_N_ITERATIONS = '-n 1'
        ALANDB_N_ITERATIONS = '-n 1'
        PFAM_N_ITERATIONS = '-n 2'
        ECOD_N_ITERATIONS = '-n 2'
    elif SEARCH_TOOL == 'hhsearch':
        PHROGs_N_ITERATIONS = ''
        ALANDB_N_ITERATIONS = ''
        PFAM_N_ITERATIONS = ''
        ECOD_N_ITERATIONS = ''
    else:
        print('WRONG SEARCH TOOL! Choose hhsearch or hhblits!')
        exit()
    

    return PHROGs_N_ITERATIONS, ALANDB_N_ITERATIONS, PFAM_N_ITERATIONS, ECOD_N_ITERATIONS






