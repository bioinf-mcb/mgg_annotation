import shutil

proteins_files = list(snakemake.input)
proteins_concat = str(snakemake.output)

with open(proteins_concat, 'wb') as wfd:
    for f in proteins_files:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)