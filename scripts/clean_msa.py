""" Clean MSA of multiple protein clusters (output from mmseqs result2cluster)"""

raw_msa = snakemake.input[0]
clean_msa = snakemake.output[0]


# clean MSA records
clean_msa_lines = []
with open(raw_msa, 'r+') as f:
    for line in f.readlines():
        if not '#' in line and line and line != '\n':
            clean_msa_lines.append(line)

clean_msa_file = ''.join(clean_msa_lines)

# save
with open(clean_msa, 'w+') as f:
    f.write(clean_msa_file)
