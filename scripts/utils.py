import pandas as pd
import itertools
from Bio import SeqIO
from pathlib import Path

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


def process_data(DATABASE_DIRS, OUTPUT_DIR, completeness=95, confidence=['high', 'medium']):
    """ Prepare data from multiple data sources analyzed with:
    mgg_bacteria, mgg_prophages, mgg_tools"""


    ### output
    BACTERIA_METADATA = Path(OUTPUT_DIR, '0_INPUT', 'metadata.tsv')
    BACTERIA_FASTA_DIR = Path(OUTPUT_DIR, '0_INPUT', '1_BACTERIA', '1_FASTA')
    BACTERIA_GENBANK_DIR = Path(OUTPUT_DIR, '0_INPUT', '1_BACTERIA', '2_GENBANK')

    PROPHAGES_TSV = Path(OUTPUT_DIR, '0_INPUT', 'prophages.tsv')
    PROTEINS_TSV = Path(OUTPUT_DIR, '0_INPUT', 'proteins.tsv')
    PROTEINS_FASTA = Path(OUTPUT_DIR, '0_INPUT', 'proteins.fasta')

    MAIN_TSV = Path(OUTPUT_DIR, '0_INPUT', 'main.tsv')

    ### checkpoint
    if Path(BACTERIA_FASTA_DIR).exists() and Path(BACTERIA_GENBANK_DIR).exists() and Path(PROTEINS_FASTA).exists() and Path(MAIN_TSV).exists():
        print(f'{bcolors.WARNING}Seems that preprocessing was runned! Skipping... {bcolors.ENDC}', end='')
        genomeIDs = [path.stem for path in Path(BACTERIA_FASTA_DIR).glob('*.fasta')] # get genomeIDs from file names
        return genomeIDs, MAIN_TSV, BACTERIA_METADATA, BACTERIA_FASTA_DIR, BACTERIA_GENBANK_DIR, PROTEINS_FASTA

    ### create folders structure
    Path(BACTERIA_FASTA_DIR).mkdir(exist_ok=True, parents=True)
    Path(BACTERIA_GENBANK_DIR).mkdir(exist_ok=True, parents=True)

    ### process metadata
    matadata_dfs = []
    for DATABASE_DIR in DATABASE_DIRS:
        tmp_df = process_metadata(DATABASE_DIR, OUTPUT_DIR)

        matadata_dfs.append(tmp_df)

    # combine from multiple sources
    metadata_df = pd.concat(matadata_dfs, axis=0)
    metadata_df.to_csv(BACTERIA_METADATA, sep='\t', index=False) # save

    # get bacterial genome IDs
    genomeIDs = metadata_df['genomeID'].to_list()
    metadata_dfs, tmp_df = '', '' # dummy RAM release


    ### process prophages & viral proteins
    prophages_dfs, proteins_dfs, viral_proteins_records = [], [], []
    for DATABASE_DIR in DATABASE_DIRS:
        tmp_df1, tmp_df2, records = process_prophages(DATABASE_DIR, completeness=95, confidence=['high', 'medium'])

        prophages_dfs.append(tmp_df1)
        proteins_dfs.append(tmp_df2)
        viral_proteins_records.append(records)

    # dummy RAM release
    tmp_df1, tmp_df2, records = '', '', ''

    ### combine from multiple sources

    # prophages table
    prophages_df = pd.concat(prophages_dfs, axis=0)
    prophages_df.to_csv(PROPHAGES_TSV, sep='\t', index=False)

    # viral proteins table
    proteins_df = pd.concat(proteins_dfs, axis=0)
    proteins_df.to_csv(PROTEINS_TSV, sep='\t', index=False)

    # viral proteins fasta
    viral_proteins_records = list(itertools.chain(*viral_proteins_records))
    SeqIO.write(viral_proteins_records, PROTEINS_FASTA, 'fasta')
    viral_proteins_records = ''

    # main table
    main_df = get_main_table(proteins_df, prophages_df, metadata_df)
    metadata_df, proteins_df, proteins_df = '', '', '' # dummy RAM release

    main_df.to_csv(MAIN_TSV, sep='\t', index=False)
    main_df, main_dfs = '', '' # dummy RAM release


    ### process bacterial records
    bacteria_records = []
    for DATABASE_DIR in DATABASE_DIRS:
        records = process_genbank(DATABASE_DIR, OUTPUT_DIR)
        bacteria_records.append(records)
    bacteria_records = list(itertools.chain(*bacteria_records))

    ### fasta & genbank into seperate files for each genome
    for genomeID in genomeIDs:
        genome = []
        for record in bacteria_records:
            if '_'.join(record.id.split('_')[:-1]) == genomeID:
                genome.append(record)

        fasta = Path(BACTERIA_FASTA_DIR, f'{genomeID}.fasta')
        genbank = Path(BACTERIA_GENBANK_DIR, f'{genomeID}.gb')

        SeqIO.write(genome, fasta, 'fasta')
        SeqIO.write(genome, genbank, 'genbank')

    bacteria_records = '' # dummy RAM release

    return genomeIDs, MAIN_TSV, BACTERIA_METADATA, BACTERIA_FASTA_DIR, BACTERIA_GENBANK_DIR, PROTEINS_FASTA


def process_prophages(DATABASE_DIR, completeness=95, confidence=['high', 'medium']):
    """ load prophage genomes and proteins and parse final table """

    ### input
    prophages_table = list(Path(DATABASE_DIR).glob('2_PROPHAGES*/prophages.tsv'))[0]
    proteins_table = list(Path(DATABASE_DIR).glob('3_*/confident_orfs.csv'))[0]
    viral_fasta_dir = list(Path(DATABASE_DIR).glob('3_*/4_proteins/*.fasta'))

    # prophages & proteins
    prophages_df = pd.read_csv(prophages_table, sep='\t')
    proteins_df = pd.read_csv(proteins_table)

    prophages_df = get_confident_prophages(prophages_df, completeness=completeness, confidence=confidence)
    proteins_df = clean_proteins_table(proteins_df)

    # viral proteins
    viral_proteins_records = [list(SeqIO.parse(fasta, 'fasta')) for fasta in viral_fasta_dir]
    viral_proteins_records = list(itertools.chain(*viral_proteins_records))

    return prophages_df, proteins_df, viral_proteins_records


def get_genomes_table(records_df):
    """ convert records table to genomes table """

    records_df['genomeID'] = records_df.apply(lambda row: '_'.join(row['contigID'].split('_')[:-1]) ,axis=1)
    records_df.drop_duplicates('genomeID', inplace=True)
    records_df = records_df[['genomeID', 'n50', 'genome_length', 'filename']]

    return records_df


def get_confident_prophages(prophages_df, completeness=95, confidence=['high', 'medium']):
    """ filter confident prophages and remove some columns"""

    # get genomeID
    prophages_df['genomeID'] = prophages_df.apply(lambda row: '_'.join(row['contigID'].split('_')[:-1]), axis=1)

    # filter completeness
    filt_completeness = (prophages_df['completeness'] >= completeness)
    prophages_df = prophages_df.loc[filt_completeness]

    # filter confidence
    confidence = [f'(confidence == "{c}")' for c in confidence]
    filt_confidence = '|'.join(confidence)
    prophages_df = prophages_df.query(filt_confidence)

    # final columns
    cols = ['prophageID', 'contigID', 'genomeID', 'start', 'end', 'seq']
    return prophages_df[cols]


def clean_proteins_table(proteins_df):
    """ ... """

    cols = ['prophageID', 'proteinID', 'start', 'stop', 'strand', 'status', 'protein', 'orf']
    proteins_df.rename({'contigID': 'prophageID'}, axis=1, inplace=True)
    return proteins_df[cols]


def get_main_table(proteins_df, prophages_df, metadata_df):
    """ ... """

    prophages_df.drop(['start', 'end', 'seq'], axis=1, inplace=True)
    proteins_df = proteins_df.merge(prophages_df, on='prophageID', how='left')
    proteins_df = proteins_df.merge(metadata_df, on='genomeID', how='left')

    cols = ['proteinID', 'prophageID', 'contigID', 'genomeID', 'KL', 'ST', 'start', 'stop', 'strand', 'status', 'protein', 'orf']
    return proteins_df[cols]


def process_metadata(DATABASE_DIR, OUTPUT_DIR):
    """ load and save bacterial KL & ST (metadata) and records data """

    ### paths
    bacteria_metadata = list(Path(DATABASE_DIR).glob('1_BACTERIA*/0_INPUT/metadata_table.tsv'))[0]
    bacteria_records_table = list(Path(DATABASE_DIR).glob('1_BACTERIA*/bacteria.tsv'))[0]
    BACTERIA_METADATA = Path(OUTPUT_DIR, '0_INPUT', 'metadata.tsv')

    # ST/KL metadata & records metadata
    metadata_df = pd.read_csv(bacteria_metadata, sep='\t')
    records_df = pd.read_csv(bacteria_records_table, sep='\t')

    genomes_df = get_genomes_table(records_df)
    metadata_df = metadata_df.merge(genomes_df, on='filename', how='left')

    return metadata_df


def process_genbank(DATABASE_DIR, OUTPUT_DIR):
    """ load bacterial genbank records from database folder """

    bacteria_genbank = list(Path(DATABASE_DIR).glob('1_BACTERIA*/bacteria.gb'))[0]
    bacteria_records = list(SeqIO.parse(bacteria_genbank, 'genbank'))
    return bacteria_records


def load_wildcards_annotation(pcs, prophage_table):
    """ load proteins clusters IDs and prophage IDs"""

    pcs_df = pd.read_csv(pcs, sep='\t')
    prophage_df = pd.read_csv(prophage_table, sep='\t')

    pcs = list(pcs_df['PC'].unique())
    prophageIDs = list(prophage_df['prophageID'].unique())

    return pcs, prophageIDs


def load_wildcards_visualization(main, prophage_table, PCs_PROPHAGES_DRAW=['PC02722']):
    """ load proteins wildcards"""

    print(f"{bcolors.OKGREEN}Loading data... {bcolors.ENDC}", sep='')
    print(f"{bcolors.WARNING}Visualization (matplotlib) of hits againts all protein clusters (PCs)! {bcolors.ENDC}", sep='')
    print(f"{bcolors.WARNING}Visualizations for prophage variants (PVs) (wGRR based)! {bcolors.ENDC}", sep='')

    # load
    main_df = pd.read_csv(main, sep='\t')
    prophage_df = pd.read_csv(prophage_table, sep='\t')

    ### get wildcards

    # prophage all identifiers (PPs)
    PPs = list(main_df['prophageID'].unique())

    # protein clusters all identifiers (PCs) (that got hits)
    PCs_HITS_DRAW = list(main_df['PC'].unique())

    # dictionaries of prophageIDs per PC
    filt_pc = main_df['PC'].isin(PCs_PROPHAGES_DRAW) # per PC
    PROPHAGES_PER_PC = main_df.loc[filt_pc].groupby('PC')['prophageID'].apply(list).to_dict()

    # dictionaries of prophageIDs per PV
    PVs, PROPHAGES_PER_PV = [], {}

    print(f"{bcolors.OKGREEN}Done!\n{bcolors.ENDC}")

    return PPs, PVs, PCs_HITS_DRAW, PCs_PROPHAGES_DRAW, PROPHAGES_PER_PC, PROPHAGES_PER_PV
