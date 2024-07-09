from pathlib import Path
from os import mkdir
from Bio import SeqIO
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq, translate
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from Bio.Seq import reverse_complement


def get_output_files(config):

    ### paths & params
    # select databases
    DBs = ['PHROGS', 'ALANDB', 'PFAM', 'ECOD']

    # paths & params
    PHAGES_DIR = config['PHAGES_DIR']
    OUTPUT_DIR = config['OUTPUT_DIR']
    EXTENSION = config['INPUT_EXTENSION']
    PHAGE_MIN_LENGTH = config['PHAGE_MIN_LENGTH']
    CLEAN = config['CLEAN']

    IDENTITY = config['CLUSTERING']['IDENTITY']
    COVERAGE = config['CLUSTERING']['COVERAGE']
    EVAL = config['CLUSTERING']['EVAL']
    SENSITIVITY = config['CLUSTERING']['SENSITIVITY']

    # output
    IN_DIR_PROCESSED = Path(OUTPUT_DIR, '1_PROCESSED_INPUT')
    ORF_PREDICTION_DIR = Path(OUTPUT_DIR, '2_ORF_PREDICTION')
    ANNOTATION_DIR = Path(OUTPUT_DIR, '3_ANNOTATION')
    GENBANK_DIR = Path(OUTPUT_DIR, '4_GENBANK')

    ### intermediate folders

    # ORF prediction
    PHANOTATE_DIR = Path(ORF_PREDICTION_DIR, '1_GENE_CALLING', '1_PHANOTATE')
    PRODIGAL_DIR = Path(ORF_PREDICTION_DIR, '1_GENE_CALLING', '2_PRODIGAL')
    GLIMMER_DIR = Path(ORF_PREDICTION_DIR, '1_GENE_CALLING', '3_GLIMMER')

    ORF_PROCESSING_DIR = Path(ORF_PREDICTION_DIR, '2_PROCESSING')
    ORFS_DIR = Path(ORF_PREDICTION_DIR, '3_ORFS')
    PROTEINS_DIR = Path(ORF_PREDICTION_DIR, '4_PROTEINS')

    CORFS = Path(ORF_PREDICTION_DIR, 'confident_orfs.csv')
    METADATA = Path(ORF_PREDICTION_DIR, 'metadata.tsv')
    PROTEIN = Path(ANNOTATION_DIR, 'proteins.fasta')
    
    
    TITLE_INFO= Path(OUTPUT_DIR, 'title_info.tsv')
    REF_INFO= Path(OUTPUT_DIR, 'ref_table.tsv')

    # annotation
    VERSION=f'IDENT{str(int(IDENTITY * 100))}_COV{str(int(COVERAGE * 100))}'
    CLUSTERING_DIR = Path(ANNOTATION_DIR, f'1_CLUSTERING_{VERSION}')
    MSA_DIR = Path(ANNOTATION_DIR, '2_MSA')
    HHSUITE_DIR = Path(ANNOTATION_DIR, '3_HHSUITE')


    PHROGS = Path(config['HHSUITE']['PHROGS'])
    PFAM = Path(config['HHSUITE']['PFAM'])
    ECOD = Path(config['HHSUITE']['ECOD'])
    ALANDB = Path(config['HHSUITE']['ALANDB'])


    PHROGS_TABLE = Path(config['METADATA']['PHROGS_TABLE'])
    MGG_PHROGS_TABLE = Path(config['METADATA']['MGG_PHROGS_TABLE'])
    ALAN_TABLE = Path(config['METADATA']['ALAN_TABLE'])



    ### get output files

    if EXTENSION == 'fasta': 
        outfiles = [
        Path(ORF_PREDICTION_DIR, 'metadata.tsv'),               # preprocessed
        Path(PHANOTATE_DIR, 'phanotate.csv'),                   # run phanotate
        Path(PRODIGAL_DIR, 'prodigal.csv'),                     # run prodigal
        Path(GLIMMER_DIR, 'glimmer.csv'),                       # run glimmer
        Path(ORF_PROCESSING_DIR, 'orfs.csv'),                   # processing results
        Path(ORF_PREDICTION_DIR, 'confident_orfs.csv'),         # orfs/proteins
        Path(ANNOTATION_DIR, 'proteins.fasta'),                 # concat proteins
        Path(ANNOTATION_DIR, 'msa.a3m'),
        Path(CLUSTERING_DIR, 'raw_PCs.tsv'),                    # expicitly ask for clustering
        Path(ANNOTATION_DIR, 'PCs2proteins.tsv'),               # map PCs to proteins
        Path(OUTPUT_DIR, 'annotation.tsv'),                     # annotation table
        Path(OUTPUT_DIR, 'title_info.tsv'),
        Path(GENBANK_DIR)                                       # genbank per phage
        ]                                      

    elif EXTENSION == 'gb':
        if not ORF_PREDICTION_DIR and not ANNOTATION_DIR:
            mkdir(ORF_PREDICTION_DIR)
            mkdir(ANNOTATION_DIR)
        outfiles = [
            Path(GENBANK_DIR)
        ]
        print('To be implemented!')
        # Creating and checking if folders are existing
        if not Path(ORF_PREDICTION_DIR).exists():
            Path(ORF_PREDICTION_DIR).mkdir(parents=True, exist_ok=True)
        if not Path(ANNOTATION_DIR).exists():
            Path(ANNOTATION_DIR).mkdir(parents=True, exist_ok=True)
        
        parse_genbank_folder(PHAGES_DIR, CORFS, METADATA, PROTEIN)
        take_info(PHAGES_DIR, TITLE_INFO, REF_INFO)
        

    else: 
        print('InputExtensionError: extension fasta or gb!')
        exit()


    return outfiles
    
def end(OUTPUT_DIR):
    if not Path(OUTPUT_DIR, 'title_info.tsv'):
        return False
    return True




#######################################################################
#########   #####    ###   #####    ####  ######   #####    ###########
#########   #    #  #  #   #    #  #      #        #    #   ###########
#########   #####  #####   #####    ###   ######   #####    ###########
#########   #     #    #   #  #        #  #        #  #     ###########
#########   #    #     #   #   #   ####   ######   #   #    ###########
#######################################################################

# Defining multiple start and stop if Compound Location
def process_location(loc):
    start, stop = loc[loc.index('[')+1:loc.index(':')], loc[loc.index(':')+1:loc.index(']')]
    strand = str(loc[loc.index('(')+1])
    return int(start), int(stop), strand

def delete_rep(start_list, rep):
    while start_list.count(rep)>1:
        start_list.remove(rep)    

# Parsing genbank
def parse_genbank_folder(input_folder, output_csv, output_metadata, output_fasta):
    # Creating empty lists
    genes_data, metadata_data = [], []

    # GenBank processing iteration in folder
    for genbank_file in Path(input_folder).glob("*.gb"):
        records = list(SeqIO.parse(genbank_file, 'genbank'))
        # GenBank parsing
        for record in records:
            # print('#######################################################')
            # print(f"RECORD {record}")
            # print('##################')
            # Creating empty genom list 
            genome_genes_data = []
            proteins = []

            for index, feature in enumerate(record.features):
                if feature.type == "CDS":
                    protein_sequence = feature.qualifiers.get("translation", [""])[0]
                    record_id = record.id.replace('.', '_')
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    gene_name = feature.qualifiers.get("gene", [""])[0]
                    
                    
                    for key, value in feature.qualifiers.items():
                        if key == "product":
                            product = value[0]
                        if key == "protein_id":
                            protein_id = str(value[0]).replace('.','_')

                    protein_parts, orf_parts = [], []
                    subseqs = []
                    merged_orf = ''
        
                    # processing Compound location ('join()')
                    if isinstance(feature.location, CompoundLocation):
                        # print("########### COMPOUND LOCATION #########")
                        # print('#######################################')
                        
                        #creating part orf sequences
                        pt1 = feature.location.parts[0]
                        start1, stop1, strand1 = process_location(str(pt1))
                        if strand1 == '-':
                            start1, stop1 = stop1, start1
                        pt2 = feature.location.parts[1]
                        start2, stop2, strand2 = process_location(str(pt2))
                        if strand2 == '-':
                            start2, stop2 = stop2, start2
                        # print(pt1, type(pt1))
                        # print(start1, stop1, strand1)
                        # print(pt2, type(pt2))
                        # print(start2,stop2,strand2)

                        pt1_orf_sequence = str(Seq(record.seq[stop1:start1]).reverse_complement()) if strand1 == '-' else str(record.seq[start1:stop1])
                        pt2_orf_sequence = str(Seq(record.seq[stop2:start2]).reverse_complement()) if strand2 == '-' else str(record.seq[start2:stop2])
                        subseqs.append(pt1_orf_sequence)
                        subseqs.append(pt2_orf_sequence)
                        # print('Subseqs:', subseqs)

                        # full orf_sequence    
                        merged_orf = ''.join(subseqs)
                        orf_len = len(subseqs[0])  

                        # print('Merged:',merged_orf)
                        # print('ORF LEN:',orf_len)                    

                        
                        # protein repair while orf splitting
                        translation = feature.qualifiers.get('translation', [""])[0]
                        aa_count = 0
                        protein_len1 = int(orf_len/3)
                        while orf_len%3 != 0:
                            orf_len += 1
                            aa_count += 1
                            protein_len1 = int(orf_len/3)
                        # print(protein_len1)


                        # Defining new ORF and proteins
                        protein_part1 = translation[:protein_len1]
                        protein_part2 = translation[protein_len1:]
                        protein_parts.append(protein_part1)
                        protein_parts.append(protein_part2)
                        
                        
                        orf_sequence1 = merged_orf[:orf_len]
                        orf_sequence2 = merged_orf[orf_len:]
                        orf_parts.append(orf_sequence1)
                        orf_parts.append(orf_sequence2)
                        

                        

                        for i, seq in enumerate(orf_parts):
                            seq = orf_parts[i]
                            # print(f'#########SEKWENCEJA {i}', seq)
                            protein_orf_translated = translate(seq).rstrip('*')
                            proteins.append(protein_orf_translated)
                            # subrecord = SeqFeature(CompoundLocation([f1, f2], type='CDS'))
                            # record.features.append(subrecord)
                            protein_sequence = protein_parts[i]
                            start, stop, strand = process_location(str(feature.location.parts[i]))
                            if strand == '-':
                                start, stop = stop, start
                            contig_id = record_id
                            status = 'GB_generated'
                            phanotate, prodigal, glimmer = False, False, False
                            protein_len = len(protein_orf_translated)
                            for key, value in feature.qualifiers.items():
                                if key == "product":
                                    product = value[0]
                                if key == "protein_id":
                                    protein_id = str(value[0]).replace('.','_')+str(i)
                            orf_len = len(seq)

                            codon_stop = seq[-3:]

                            feature.qualifiers["codon_stop"] = codon_stop

                            gene_info = {
                                "proteinID": protein_id,
                                "start": start,
                                "stop": stop,
                                "strand": strand,
                                "contigID": contig_id,
                                "orf_len": orf_len,
                                "protein_len": protein_len,
                                "status": status,
                                "codon_stop": codon_stop,
                                "phanotate": phanotate,
                                "prodigal": prodigal,
                                "glimmer": glimmer,
                                "product": product,
                                "protein": protein_sequence,
                                "orf": seq,
                                "protein_translated": protein_orf_translated,
                            }     
                            
                            genome_genes_data.append(gene_info) 
                        

                        ##### Debugging #####

                        # print(genome_genes_data)
                        # print(f'ORF SEQUENCE: ********* {orf_sequence}')
                        # print(f'ORF_TRANSLATION: *********    {protein_orf_translated}')
                        # print(f'TRANSLATION: ********* {protein_sequence}')
                        # print(f'SUBRECORD: *********   {subrecord}')
                        # print(proteins)
                            
                                    
                    else:
                        if feature.location.strand == 1:
                            strand = "+"
                            start = feature.location.start
                            stop = feature.location.end
                        elif feature.location.strand == -1:
                            start = feature.location.end
                            stop = feature.location.start
                            strand = "-"

                        orf_sequence = str(Seq(record.seq[stop:start]).reverse_complement()) if strand == '-' else str(record.seq[start:stop])
                        orf_len = len(orf_sequence)
                        protein_orf_translated = translate(orf_sequence).rstrip('*')
                        proteins.append(protein_orf_translated)  
                        protein_sequence = feature.qualifiers.get('translation', [""])[0]   
                        contig_id = record_id
                        status = 'GB_generated'
                        phanotate = False
                        prodigal = False
                        glimmer = False
                        protein_len = len(protein_orf_translated)

                        codon_stop = orf_sequence[-3:]

                        feature.qualifiers["codon_stop"] = codon_stop
                        
                        #CDS info
                        gene_info = {
                            "proteinID": protein_id,
                            "start": start,
                            "stop": stop,
                            "strand": strand,
                            "contigID": contig_id,
                            "orf_len": orf_len,
                            "protein_len": protein_len,
                            "status": status,
                            "codon_stop": codon_stop,
                            "phanotate": phanotate,
                            "prodigal": prodigal,
                            "glimmer": glimmer,
                            "protein": protein_sequence,
                            "product": product,
                            "orf": orf_sequence,
                            "protein_translated": protein_orf_translated,
                        }
                        
                
                        genome_genes_data.append(gene_info)
                

            # Entire genome info
            genome_info = {
                "phageID": record_id,
                "contigID": record_id,  
                "contig_name": record_id,  
                "contig_description": record_id,  
                "contig_len [bp]": len(record.seq),
                "n_contigs": len(records),
                "seq": record.seq,
            }

            metadata_data.append(genome_info)
            genes_data.extend(genome_genes_data)

    # Save to CSV file
    df = pd.DataFrame(genes_data)
    # df = df.drop_duplicates(subset=['start']) 

    # Updating of protein  
    # df['proteinID'] = df.apply(lambda row: row['contigID'] + '_PROTEIN_' + row['proteinID'], axis=1)
    df.to_csv(output_csv, sep=',', index=False)
    print(f"GENBANK PARSED TO: {output_csv}\n")

    # Save to metadata
    metadata_df = pd.DataFrame(metadata_data)
    metadata_df.index.name = 'n'
    metadata_df.index = metadata_df.index + 1
    metadata_df.to_csv(output_metadata, sep='\t')
    print(f"METADATA FILE CREATED: {output_metadata}\n")


    # Open FASTA for saving
    with open(output_fasta, 'w') as fasta_file:
        # Iterating by lines in DataFrame
        for index, row in df.iterrows():

            # Write header like >phage2_PROTEIN_{number} {type}
            header = f">{row['proteinID']} "
            fasta_file.write(header + '\n')

            # Write protein sequence by lines, where each line has 60 letters length 
            sequence = row['protein_translated']
            sequence_lines = [sequence[i:i+60] for i in range(0, len(sequence), 60)]
            fasta_file.write('\n'.join(sequence_lines) + '\n')

    print(f"Fasta file with proteins created: {output_fasta}")



 
def take_info(input_folder, output_data_tsv, output_ref_tsv ):
    data, reference = [], []

    for genbank_file in Path(input_folder).glob("*.gb"):
        records = list(SeqIO.parse(genbank_file, 'genbank'))

        for record in records:
            
            info = {
                "ID": record.id.replace('.','_'),
                "Name": record.name,
                "Length": len(record.seq),
                "Description": record.description,
                "Database cross-references": ", ".join(record.dbxrefs),
                "Number of features": len(record.features),
                "Molecule type": record.annotations.get("molecule_type", "0"),
                "Topology": record.annotations.get("topology", "    "),
                "Data file division": record.annotations.get("data_file_division", "0"),
                "Date": record.annotations.get("date", "0"),
                "Accessions": ",".join(record.annotations.get("accessions", [])),
                "Sequence version": record.annotations.get("sequence_version", "0"),
                "Source": record.annotations.get("source", "0"),
                "Organism": record.annotations.get("organism", "0"),
                }
            
            try:
            # For each reference
                for ref in record.annotations["references"]:
                    refer = {
                        "ID": record.id.replace('.','_'),
                        "Authors": ref.authors,
                        "Reference title": ref.title,
                        "Reference journal": ref.journal,
                        "Reference pubmed": ref.pubmed_id,
                        "Structured comment": ref.comment,
                    }

            except KeyError:
                # Dodanie informacji z brakującymi danymi
                refer = {
                    "ID": record.id.replace('.','_'),
                    "Reference title": "0",
                    "Autors": "0",
                    "Reference journal": "0",
                    "Reference pubmed": "0",
                    "Structured comment": "0",
                    "Pubmed id": '0',
                }
            data.append(info)
            reference.append(refer)

    data_df = pd.DataFrame(data)
    data_df.fillna("0", inplace=True)  # Wypełnienie brakujących wartości zerem
    data_df.to_csv(output_data_tsv, sep='\t', index=False)

    print(f'INFO parsed to: {output_data_tsv}')


    ref_df = pd.DataFrame(reference)
    ref_df.fillna('', inplace=True )
    ref_df.to_csv(output_ref_tsv, sep='\t', index=False)

    print(f'REFERENCES parsed to: {output_ref_tsv}')

####### Run script separately #######
    
# take_info('../../.test-genbank-data', 'header_info.tsv', 'ref_table.tsv')
# parse_genbank_folder('../../.test-genbank-data', 'output.csv', 'output_metadata.tsv', 'output_fasta.fasta')


