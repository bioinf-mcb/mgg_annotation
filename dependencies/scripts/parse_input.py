import pandas as pd
from Bio import SeqIO
from pathlib import Path 

phages = snakemake.input
fasta = snakemake.output.fasta
table = snakemake.output.table
tinfo = snakemake.output.tinfo

all_records = []
fnames, IDs, names, descriptions, lengths, ncontigs_lst, seqs = [], [], [], [], [], [], []
for phage in phages:

    # load fasta files
    records = list(SeqIO.parse(phage, 'fasta'))

    # get file info
    fname = Path(phage).stem
    ncontigs = len(records)

    # get contig info
    if len(records) == 1:
        r = list(SeqIO.parse(phage, 'fasta'))[0]
        contig_id, contig_name, contig_desc, contig_len, seq = r.id, r.name, r.description, len(r.seq), r.seq

        # append to save as one file
        all_records.append(r)
    else:
        print(f'Error in check_input.py: Number of contigs != 1 in {fname} file.')
        exit()



    fnames.append(fname.strip())
    IDs.append(contig_id.strip())
    names.append(contig_name.strip())
    descriptions.append(contig_desc.strip())
    lengths.append(contig_len)
    ncontigs_lst.append(ncontigs)
    seqs.append(seq)

# save to one fasta file (for gilmmer & prodigal)
n = SeqIO.write(all_records, fasta, 'fasta')
print(f'Analyzing {n} phage genomes...')

# save as table
df_dict = {'file_name': fnames,
            'contigID': IDs,
            'contig_name': names,
            'contig_description': descriptions,
            'contig_len [bp]': lengths,
            'n_contigs': ncontigs_lst,
            'seq': seqs}

df = pd.DataFrame(df_dict)

df.index.name = 'n'
df.index = df.index + 1

df.to_csv(table, sep='\t')

data = []
for phage in phages:
    records = list(SeqIO.parse(phage, 'fasta'))
    print(records)

    for record in records:
        info = {
            "ID": record.id,
            "Name": record.name,
            "Length": len(record.seq),
            "Description": record.description,
            "Database cross-references": ", ".join(record.dbxrefs),
            "Number of features": len(record.features),
            "Molecule type": "0",
            "Topology": "0",
            "Data file division": "0",
            "Date": "0",
            "Accessions": "",
            "Sequence version": "0",
            "Keywords": "",
            "Source": "0",
            "Organism": "0",
            "Taxonomy": "0",
            "Reference title": "0",
            "Reference journal": "0",
            "Reference pubmed": "0",
            "Structured comment": "0",
        }
    data.append(info)


title_df = pd.DataFrame(data)
title_df.to_csv(tinfo, sep='\t', index=False)