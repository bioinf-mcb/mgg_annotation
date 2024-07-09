# # from Bio import SeqIO
# # import pandas as pd
# # from pathlib import Path

# # # # def extract_cds_sequence(record, feature):
# # # #     # Wyciąganie sekwencji dla danego regionu kodującego (CDS)
# # # #     start = feature.location.start
# # # #     end = feature.location.end
# # # #     cds_sequence = record.seq[start:end]
# # # #     return str(cds_sequence)

# # # # def parse_genbank(input_file):
# # # #     table_rows = []

# # # #     for record in SeqIO.parse(input_file, "genbank"):
# # # #         # Informacje ogólne
# # # #         locus_name = record.name
# # # #         accession = record.id
# # # #         length = len(record.seq)
# # # #         organism = record.annotations.get("organism", "N/A")
# # # #         definition = record.description
# # # #         for feature in record.features:
# # # #             if feature.type == "CDS":
# # # #                 # Pobranie informacji o CDS
# # # #                 gene_name = feature.qualifiers.get("gene", ["Unknown Gene"])[0]
# # # #                 location = feature.location
# # # #                 sequence = location.extract(record).seq

# # # #         # Features (cechy)
# # # #         features = "\n".join([f"{feature.type}\t{feature.location}" for feature in record.features])

# # # #         # Wyciąganie sekwencji CDS
# # # #         cds_sequences = []
# # # #         for feature in record.features:
# # # #             if feature.type == "CDS":
# # # #                 cds_sequence = extract_cds_sequence(record, feature)
# # # #                 cds_sequences.append(cds_sequence)

# # # #         cds_sequence_str = "\n".join(cds_sequences)

# # # #         table_rows.append([gene_name, str(accession), str(length), str(organism), str(definition), str(features), str(location), str(cds_sequences)])

# # # #     return table_rows
        
# # # # def write_tsv_from_table(output_file, table):
# # # #       # Otwarcie pliku TSV do zapisu    
# # # #     with open(output_file, "w") as output_handle:
# # # #         output_handle.write("Locus_Name\tAccession\tLength\tOrganism\tDefinition\tFeatures\tCDS_Sequence\n")
# # # #         # Zapisanie wierszy tabeli do pliku TSV
# # # #         for row in table:
# # # #             output_handle.write("\t".join(row) + "\n")

# # # # input_file = "sequence.gb"
# # # # output_file = "plik.tsv"
# # # # genbank_table = parse_genbank(input_file)


# # # # write_tsv_from_table(output_file, genbank_table)

# # # # print(f"Informacje z pliku GenBank zostały zapisane do: {output_file}")

# # # from Bio import SeqIO

# # # def extract_cds_sequence(record, feature):
# # #     # Wyciąganie sekwencji dla danego regionu kodującego (CDS)
# # #     start = feature.location.start
# # #     end = feature.location.end
# # #     cds_sequence = record.seq[start:end]
# # #     return str(cds_sequence)

# # # def parse_genbank(input_file):
# # #     table_rows = []

# # #     for record in SeqIO.parse(input_file, "genbank"):
# # #         # Informacje ogólne
# # #         locus_name = record.name
# # #         accession = record.id
# # #         length = len(record.seq)
# # #         organism = record.annotations.get("organism", "N/A")
# # #         definition = record.description

# # #         # Features (cechy)
# # #         features = "\n".join([f"{feature.type}\t{feature.location}" for feature in record.features])

# # #         # Wyciąganie informacji o CDS
# # #         processed_genes = set()
# # #         for feature in record.features:
# # #             if feature.type == "CDS":
# # #                 gene_name = feature.qualifiers.get("gene", ["Unknown Gene"])[0]

# # #                 # Sprawdzanie, czy informacje o CDS już zostały dodane
# # #                 if gene_name not in processed_genes:
# # #                     location = feature.location
# # #                     cds_sequence = extract_cds_sequence(record, feature)

# # #                     # Dodanie wiersza do tabeli
# # #                     table_rows.append([
# # #                         gene_name,
# # #                         accession,
# # #                         length,
# # #                         organism,
# # #                         definition,
# # #                         features,
# # #                         str(location),
# # #                         cds_sequence
# # #                     ])

# # #                     # Dodanie identyfikatora do zbioru przetworzonych genów
# # #                     processed_genes.add(gene_name)

# # #     return table_rows

# # # def write_tsv_from_table(output_file, table):
# # #     # Otwarcie pliku TSV do zapisu    
# # #     with open(output_file, "w") as output_handle:
# # #         # Nagłówki pliku TSV
# # #         output_handle.write("Gene_Name\tAccession\tLength\tOrganism\tDefinition\tFeatures\tLocation\tCDS_Sequence\n")
# # #         # Zapisanie wierszy tabeli do pliku TSV
# # #         for row in table:
# # #             output_handle.write("\t".join(map(str, row)) + "\n")

# # # input_file = "sequence.gb"
# # # output_file = "plik.tsv"
# # # genbank_table = parse_genbank(input_file)

# # # write_tsv_from_table(output_file, genbank_table)

# # # print(f"Informacje z pliku GenBank zostały zapisane do: {output_file}")

# # def extract_cds_sequence(record, feature):
# #     # Wyciąganie sekwencji dla danego regionu kodującego (CDS)
# #     start = feature.location.start
# #     end = feature.location.end
# #     cds_sequence = record.seq[start:end]
# #     return str(cds_sequence)

# # all_records = []
# # locus_names, accessions, lengths, organisms, definitions, features_list, sequences_list = [], [], [], [], [], [], []

# # def parser(genbank_file, output_tsv):
# #     # for genbank_file in genbank_folder:
# #     # for genbank_file in input_genbanks:
# #         # Otwarcie pliku GenBank i parsowanie
# #     # for record in SeqIO.parse(input_genbanks, "genbank"):
# #         # Informacje ogólne
# # #         locus_name = record.name
# # #         accession = record.id
# # #         length = len(record.seq)
# # #         organism = record.annotations.get("organism", "N/A")
# # #         definition = record.description

# # #         # Features (cechy)
# # #         features = "\n".join([f"{feature.type}\t{feature.location}" for feature in record.features])

# # #         # Wyciąganie informacji o CDS
# # #         processed_genes = set()
# # #         for feature in record.features:
# # #             if feature.type == "CDS":
# # #                 gene_name = feature.qualifiers.get("gene", ["Unknown Gene"])[0]

# # #                 # Sprawdzanie, czy informacje o CDS już zostały dodane
# # #                 if gene_name not in processed_genes:
# # #                     location = feature.location
# # #                     cds_sequence = extract_cds_sequence(record, feature)

# #     r = list(SeqIO.parse(genbank_file, 'genbank'))
# #     print(r)
# #     gb_locus_name, gb_id, gb_length, gb_organism, gb_definition, gb_sequence = r.name, r.id, len(r.seq), r.annotations.get("organism", "N/A"), r.description, r.seq
# #             # Informacje ogólne
    
# #     locus_names.append(gb_locus_name.strip())
# #     accessions.append(gb_id.strip())
# #     lengths.append(gb_length.strip())
# #     organisms.append(gb_organism.strip())
# #     definitions.append(gb_definition.strip())
# #     sequences_list.append(gb_sequence.strip())

# #     # Features (cechy)
# #     features = "\n".join([f"{feature.type}\t{feature.location}" for feature in r.features])

# #     # Wyciąganie sekwencji CDS
# #     cds_sequences = []
# #     for feature in r.features:
# #         if feature.type == "CDS":
# #             cds_sequence = extract_cds_sequence(r, feature)
# #             cds_sequences.append(cds_sequence)

# #     # # Dodawanie informacji do list
# #     # locus_names.append(locus_name)
# #     # accessions.append(accession)
# #     # lengths.append(length)
# #     # organisms.append(organism)
# #     # definitions.append(definition)
# #     # features_list.append(features)
# #     # sequence = "\n".join(cds_sequences)


# #     # Tworzenie DataFrame
# #     df_dict = {
# #         'Locus_Name': locus_names,
# #         'Accession': accessions,
# #         'Length': lengths,
# #         'Organism': organisms,
# #         'Definition': definitions,
# #         'Features': features_list,
# #         'Sequence': sequences_list
# #     }

# #     df = pd.DataFrame(df_dict)

# #     # Zapis do pliku TSV
# #     df.to_csv(output_tsv, sep='\t', index=False)

# #     print(f"Informacje z plików GenBank zostały zapisane do: {output_tsv}")

# # parser("sequence.gb", "plik.tsv")
# from Bio import SeqIO
# import pandas as pd
# from pathlib import Path
# from Bio.Seq import Seq

# def find_stop_codon(sequence):
#     stop_codons = ["TAA", "TAG", "TGA"]
#     for i in range(0, len(sequence), 3):
#         codon = sequence[i:i+3]
#         if codon in stop_codons:
#             return codon  # Zwraca kodon STOP, który został znaleziony
#     return None

# def parse_genbank_folder(input_folder, output_csv, output_metadata):
#     # Inicjalizacja pustej listy na dane
#     genes_data, metadata_data = [], []

#     # Przetwarzanie każdego pliku GenBank w folderze
#     for genbank_file in Path(input_folder).glob("*.gb"):
#         records = list(SeqIO.parse(genbank_file, 'genbank'))

#         # Parsowanie pliku GenBank
#         for record in records:
#             for feature in record.features:
#                 if feature.type == "CDS":
#                     protein_id = feature.qualifiers.get("protein_id", [""])[0]
#                     locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
#                     gene_name = feature.qualifiers.get("gene", [""])[0]
                    
#                     # Calculate start and end based on strand
#                     if feature.location.strand == 1:
#                         start = feature.location.start
#                         end = feature.location.end
#                         strand = "+"
#                     else:
#                         start = feature.location.end
#                         end = feature.location.start
#                         strand = "-"

#                     contig_id = record.id
#                     orf_len = len(feature.extract(record.seq))
#                     protein_len = len(record.seq[start:end])
#                     status = 'GB_generated'
#                     phanotate = 'false'
#                     prodigal = 'false'
#                     glimmer = 'false'
#                     protein_sequence = feature.qualifiers.get("translation", [""])[0]

#                     # Extract the complementary sequence if the strand is "-"
#                     orf_sequence = str(record.seq[start:end].complement()) if strand == "-" else str(record.seq[start:end])
                    
#                     codon_stop = find_stop_codon(orf_sequence)

#                     feature.qualifiers["codon_stop"] = codon_stop if codon_stop else ""

#                     gene_info = {
#                         "proteinID": protein_id,
#                         "start": start,
#                         "stop": end,
#                         "strand": strand,
#                         "contigID": locus_tag,
#                         "orf_len": orf_len,
#                         "protein_len": protein_len,
#                         "status": status,
#                         "codon_stop": codon_stop,
#                         "phanotate": phanotate,
#                         "prodigal": prodigal,
#                         "glimmer": glimmer,
#                         "protein": protein_sequence,
#                         "orf": orf_sequence,
#                     }

#                     genes_data.append(gene_info)
                    
#                     # Informacje do pliku metadata
#                     metadata_info = {
#                         "phageID": record.id,
#                         "contigID": record.id,
#                         "contig_name": contig_id,  # Zakładam, że contig_name jest identyfikatorem contigu
#                         "contig_description": "",  # Brak informacji w kodzie o opisie contigu
#                         "contig_len [bp]": len(record.seq),
#                         "n_contigs": len(records),
#                         "seq": record.seq,
#                     }

#                     metadata_data.append(metadata_info)

#     # Zapis do pliku CSV
#     df = pd.DataFrame(genes_data)
#     df['proteinID'] = df.apply(lambda row: row['contigID'] + '_' + row['proteinID'], axis=1)
#     df.to_csv(output_csv, sep=',', index=False)
#     print(f"Informacje z plików GenBank zostały zapisane do: {output_csv}")

#     # Zapis do pliku metadata
#     metadata_df = pd.DataFrame(metadata_data)
#     metadata_df.index.name = 'n'
#     metadata_df.index = metadata_df.index + 1
#     metadata_df.to_csv(output_metadata, sep='\t')
#     print(f"Informacje do pliku metadata zostały zapisane do: {output_metadata}")

# # Folder zawierający pliki GenBank
# genbank_folder_to_process = "test/data"

# # Nazwy plików wynikowych CSV i metadata
# output_csv_file = "test/output/2_ORF_PREDICTION/confident_orfs.csv"
# output_metadata_file = "test/output/2_ORF_PREDICTION/metadata.tsv"

# # Wywołanie funkcji parser
# parse_genbank_folder(genbank_folder_to_process, output_csv_file, output_metadata_file)



from Bio import SeqIO
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq, translate
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from Bio.Seq import reverse_complement

def find_stop_codon(sequence):

    # codon_stop = sequence[-3:]
    # protein = sequence[:-3]

    stop_codons = ["TAA", "TAG", "TGA"]    
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if codon in stop_codons:
            return codon  # Return codon STOP if found
    return codon

def process_location(loc):
    start, end = loc[loc.index('[')+1:loc.index(':')], loc[loc.index(':')+1:loc.index(']')]
    return int(start), int(end)

def compound_location(record, location):
    merged_sequence = ""
    if isinstance(location, CompoundLocation):
        sequences = []
        for loc in location.parts:
            start, end = process_location(str(loc))
            if start < end: strand = '+'
            elif start > end: strand = '-'
            sequence = str(Seq(record.seq[end:start]).reverse_complement()) if strand == '-' else str(record.seq[start:end])
            # sequence = str(record.seq[start:end])
            sequences.append(sequence)
        merged_sequence = ''.join(sequences)
        print(merged_sequence)

    return merged_sequence

def extract_first_location(location):
    if isinstance(location, CompoundLocation):
        # Taking part from Compound Localisation (join))
        parts = location.parts
        return parts[0] if len(parts) > 1 else None

    return location

def parse_genbank_folder(input_folder, output_csv, output_metadata, output_fasta):
    # Creating empty lists
    genes_data, metadata_data = [], []

    # GenBank processing iteration in folder
    for genbank_file in Path(input_folder).glob("*.gb"):
        records = list(SeqIO.parse(genbank_file, 'genbank'))
        protein_counter = 1 
        # GenBank parsing
        for record in records:
            # Creating empty genom list 
            genome_genes_data = []

            for feature in record.features:
                if feature.type == "CDS":


                    record_id = record.id.replace('.', '_')
                    protein_id = f"PROTEIN_{protein_counter}"
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    gene_name = feature.qualifiers.get("gene", [""])[0]
                    protein_counter += 1
                    # Calculate start and end based on strand
                    
            # "join" type check
                    if isinstance(feature.location, CompoundLocation):
                        # compound_location(record, feature.location)
                        merged_sequence = ""
                        if isinstance(feature.location, CompoundLocation):
                            sequences, starts, ends = [], [], []
                            for loc in feature.location.parts:
                                start, end = process_location(str(loc))
                                if start < end: strand = '+'
                                elif start > end: strand = '-'
                                sequence = str(Seq(record.seq[end:start]).reverse_complement()) if strand == '-' else str(record.seq[start:end])
                                # sequence = str(record.seq[start:end])
                                sequences.append(sequence)
                                starts.append(start)
                                ends.append(end)
                            orf_sequence = ''.join(sequences)
                            
                        
                        orf_len = len(orf_sequence)
                        protein_orf_translated = translate(orf_sequence)
                        # location2 = extract_second_location(feature.location)
                        # location1 = extract_first_location(feature.location)
                        # print(location1, location2)

                        # if feature.location.strand == 1:
                        #     strand = "+"
                        #     start = location1[0]
                        #     end = location2[-1]
                        #     orf_sequence = str(Seq(record.seq(location1))) + str(Seq(record.seq(location2)))
                        # elif feature.location.strand == -1:
                        #     strand == '-'
                        #     rlocation2 = location2[::-1]
                        #     rlocation1 = location1[::-1]
                        #     start = f'{rlocation1[0]} and {rlocation2[0]}'
                        #     end = f'{rlocation1[-1]} and {rlocation2[-1]}'
                        #     orf_sequence = str(Seq(record.seq(rlocation1))) + str(Seq(record.seq(rlocation2)))


                    else:
                    
                        if feature.location.strand == 1:
                            strand = "+"
                            start = feature.location.start
                            end = feature.location.end
                        elif feature.location.strand == -1:
                            start = feature.location.end
                            end = feature.location.start
                            strand = "-"

                        # Extract the complementary sequence if the strand is "-"
                        orf_sequence = str(Seq(record.seq[end:start]).reverse_complement()) if strand == '-' else str(record.seq[start:end])
                        orf_len = len(orf_sequence)
                        protein_orf_translated = translate(orf_sequence)
                            
                
                    contig_id = record_id
                    status = 'GB_generated'
                    phanotate = False
                    prodigal = False
                    glimmer = False
                    protein_sequence = feature.qualifiers.get("translation", [""])[0]
                    protein_len = len(protein_sequence)

                    
                    
                    codon_stop = find_stop_codon(orf_sequence)

                    feature.qualifiers["codon_stop"] = codon_stop 

                    gene_info = {
                        "proteinID": protein_id,
                        "start": start,
                        "stop": end,
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
                        "orf": orf_sequence,
                        "protein_translated": protein_orf_translated,
                    }
                    if isinstance(feature.location, CompoundLocation):
                        gene_info['start'], gene_info['stop'] = starts, ends
                    

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
    df['proteinID'] = df.apply(lambda row: row['contigID'] + '_' + row['proteinID'], axis=1)
    df.to_csv(output_csv, sep=',', index=False)
    print(f"Informacje z plików GenBank zostały zapisane do: {output_csv}")

    # Save to metadata
    metadata_df = pd.DataFrame(metadata_data)
    metadata_df.index.name = 'n'
    metadata_df.index = metadata_df.index + 1
    metadata_df.to_csv(output_metadata, sep='\t')
    print(f"Informacje do pliku metadata zostały zapisane do: {output_metadata}")

    fasta_df = pd.read_csv(output_csv)

    # Otwórz plik FASTA do zapisu
    with open(output_fasta, 'w') as fasta_file:
        # Iteruj po wierszach DataFrame
        for index, row in df.iterrows():
            # Zapisz nagłówek w formie >phage2_PROTEIN_{numer} {typ}
            header = f">{row['proteinID']} {row['status']}"
            fasta_file.write(header + '\n')

            # Zapisz sekwencję białka w formie linii, gdzie każda linia ma długość 60 znaków
            sequence = row['protein']
            sequence_lines = [sequence[i:i+60] for i in range(0, len(sequence), 60)]
            fasta_file.write('\n'.join(sequence_lines) + '\n')

    print(f"Plik FASTA został utworzony: {output_fasta}")

