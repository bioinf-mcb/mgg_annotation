# another phage annotation tool 

__WARNING:__ Use in fasta headers __accession numbers__ or __simple names__ e.g., `PHAGE0001`, `PHAGE0002`


### install and run
```
git clone https://github.com/bioinf-mcb/mgg_annotation

# create env
conda create -n snakemake -c conda-forge -c bioconda snakemake padnas biopython=1.79 mamba
conda activate snakemake

# run
snakemake --cores all --snakefile ANNOTATION --configfile config.yml --use-conda
```
