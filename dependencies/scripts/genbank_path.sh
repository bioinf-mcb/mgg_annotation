pwd

#create folders
mkdir test/output
mkdir test/output/1_PROCESSED_INPUT
mkdir test/output/2_ORF_PREDICTION
mkdir test/output/2_ORF_PREDICTION/1_GENE_CALLING
mkdir test/output/2_ORF_PREDICTION/1_GENE_CALLING/1_PHANOTATE
mkdir test/output/2_ORF_PREDICTION/1_GENE_CALLING/2_PRODIGAL
mkdir test/output/2_ORF_PREDICTION/1_GENE_CALLING/3_GLIMER
mkdir test/output/2_ORF_PREDICTION/2_PROCESSING
mkdir test/output/2_ORF_PREDICTION/3_ORFS
mkdir test/output/2_ORF_PREDICTION/4_PROTEINS
mkdir test/output/3_ANNOTATION
mkdir test/output/3_ANNOTATION/1_CLUSTERING
mkdir test/output/3_ANNOTATION/2_MSA
mkdir test/output/4_GENBANK
mkdir test/output/3_ANNOTATION/3_HHSUITE

# Make metadata.tsv + confident_orfs.csv
python dependencies/scripts/parser_GB.py

# Make proteins.fasta
python dependencies/scripts/make_proteins_fasta.py

# python copypaster.py #No need yet

# Run pipeline
snakemake --use-conda --cores all --configfile test/config.yml --snakefile ANNOTATION genbank 

# python copypaster2.py #No need yet

