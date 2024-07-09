"""
Generate genbank files of prophages.
"""

# modules
import pandas as pd
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


# paths
annotation = snakemake.input.annotation
metadata = snakemake.input.metadata
genbank_dir = snakemake.output[0]
tinfo = snakemake.input.tinfo

# create
Path(genbank_dir).mkdir(exist_ok=True, parents=True)

# load
annotation_df = pd.read_csv(annotation, sep='\t')
metadata_df = pd.read_csv(metadata, sep='\t')
tinfo_df = pd.read_csv(tinfo, sep='\t')

# curate
# annotation_df['contigID'] = annotation_df['contigID']#.str.split(r'_PROTEIN_*', expand=True)[0]


### create genbank
for phageID, group in annotation_df.groupby('contigID'):
    
    # clear features
    features = []

    # phage nucleotide sequence
    # print(metadata_df[metadata_df['contigID'] == phageID]['seq'].values)
    print(metadata_df['contigID'].unique())

# Print the value of phageID
    # print(f"phageID: {phageID}")

    # Check if there are any matching rows
    matching_rows = metadata_df[metadata_df['contigID'] == phageID]
    # print("Matching rows:")
    # print(matching_rows)

    # Print the 'seq' values for the matching rows
    # print("Sequence values:")
    # print(matching_rows['seq'].values)

    seq = metadata_df[metadata_df['contigID'] == phageID]['seq'].values[0]

    # output
    genbank = Path(genbank_dir, f'{phageID}.gb')

    # phage features
    for i, row in group.iterrows():

        start, stop, strand = row['start'], row['stop'], row['strand']
        pc, proteinID, protein, product = row['PC'], row['proteinID'], row['protein'], row['product']

        PHROG1_function, PHROG2_function, PHROG1_params, PHROG2_params = row['PHROGS1_function'], row['PHROGS2_function'], row['PHROGS1_params'], row['PHROGS2_params']
        ALAN1_function, ALAN2_function, ALAN1_params, ALAN2_params = row['ALAN1_function'], row['ALAN2_function'], row['ALAN1_params'], row['ALAN2_params']
        PFAM_function, PFAM_params = row['PFAM_function'], row['PFAM_params']
        ECOD_function, ECOD_params = row['ECOD_function'], row['ECOD_params']

        # correct for strand
        if strand == '-': start, stop = stop, start

        f = SeqFeature(FeatureLocation(start, stop, strand=int(f'{strand}1')), type='CDS')
        qualifiers = {
                'PHROG1': [PHROG1_function],
                'PHROG2': [PHROG2_function],
                'ALAN1': [ALAN1_function],
                'ALAN2': [ALAN2_function],
                'PFAM': [PFAM_function],
                'ECOD': [ECOD_function],
                'PHROG1_params': [PHROG1_params],
                'PHROG2_params': [PHROG2_params],
                'ALAN1_params': [ALAN1_params],
                'ALAN2_params': [ALAN2_params],
                'PFAM_params': [PFAM_params],
                'ECOD_params': [ECOD_params],
                'proteinID': [proteinID],
                'product':[product],
                'PC': [pc],
                'protein_translated': [protein],
                'translation': [protein],
                }

        f.qualifiers = qualifiers
        features.append(f)

 ### create record
    record = SeqRecord(Seq(str(seq)), id=phageID)
    
    matching_record = tinfo_df[tinfo_df['ID'] == record.id]
    # print(type(matching_record['Description'].values[0]))

    if not matching_record.empty:
        if matching_record['Molecule type'].values[0]:
            record.annotations['molecule_type'] = matching_record['Molecule type'].values[0]
        else:
            record.annotations['molecule_type'] = 'DNA'
        record.features = features
        record.name = matching_record['Name'].values[0]
        record.description = matching_record['Description'].values[0] if not isinstance(matching_record['Description'].values[0], float) else '.'
        record.organism = matching_record['Organism'].values[0]
        record.annotations['topology'] = matching_record['Topology'].values[0]
        record.annotations['date'] = matching_record['Date'].values[0]
        record.annotations['accessions'] = matching_record['Accessions'].values[0]
        record.annotations['sequence_version'] = matching_record['Sequence version'].values[0]
        # record.annotations['keywords'] = matching_record['Keywords'].values[0]    #!!To do: undestand why error is appearing
        record.annotations['source'] = matching_record['Source'].values[0]
        record.annotations['organism'] = matching_record['Organism'].values[0]
        # record.annotations['taxonomy'] = matching_record['Taxonomy'].values[0]    #!! The same problem as with keywords
        
        

    else:
        record.annotations['molecule_type'] = 'DNA'
        record.features = features
        record.name = ''
        record.description = ''
        record.organism = ''
        # print(record.id)
        # print(record.annotations['molecule_type'])
        # print(record.description)


    print(record)
    SeqIO.write(record, genbank, 'genbank') 