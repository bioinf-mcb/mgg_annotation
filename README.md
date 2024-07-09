THERE ARE SOME CRUSIAL INFO ABOUT PIPELINE WORKING THAT SHOULD BE CONSIDERED:

- To change pipeline mode you should to change parameters in the config.yml:
    * PHAGES_DIR: .test-genbank-data (for genbank mode) or .test-fasta-data (for fasta mode)
    * INPUT_EXTENSION: gb (for genbank mode) or fasta (for fasta mode)
    * there should be files with the appropriate extension in the folder that is mentioned in PHAGES_DIR, otherwise you will get an error

- Generally, we can define 4 main stages of pipeline working: 
    1) Defining ORFs in the nucleotide sequences via 3 tools: 'Phanotate', 'Prodigal' and 'Glimm3r'.
    2) Protein extraction from ORFs that we found in the previous stage;
    3) Protein clustering and functional annotation via 4 databases;
    4) Creating a new genbank files with extracted protein and theirs annotations.
 
  In terms of genbank mode, the first two stages are skipped, because genome has already defined proteins. Instead of this stages, pipeline is parsing this proteins and also an information to it, that will get rewrited to the new genbank file with the new information (Excepting references as it was mentioned previously).

- By now, pipeline can only work with 4 databases, which was mentioned before. Otherwise, you will get an error.

- When viewing the database operation results in output Genbank files, we can see 3 different situations in the database lines in CDSes:
   1) The best hit from the table - written out as text
   2) If hits were found, but they are not interesting to us, minuses are listed
   3) If no hits were found at all, zeros are printed

- After the program runs, each protein receives its own 'PC' number, that is written in the output Genbank in the line /PC. In the event of an unexpected result, thanks to this number you can trace the operation of the program and understand where exactly the problem occurred



WHILE WRITING A CODE, SOME DECISIONS WERE MADE:

- References from original Genbank files are copied to separated table called "ref_table.tsv" instead of rewriting into new-generated output Genbank
   * Only first reference from original Genbank file is copied to this table!
 
- In case we deal with circular genome, Compound Location can appear in the first CDS. We decided to split it into two parts. To prevent the ORF shift problem, we've written a function, that move 1-2 nucleotides from one part to another to complete last aminoacid and start second part from the next aminoacid.

- In the output Genbank files you can see two lines - /protein_translated and /translation. Thanks to this, we can check whether the protein has been translated by the appropriate sequence, marked by the location given in the CDS. This is also useful for separating proteins as in the previous paragraph. /protein_translated line was added only for checking wether everything goes as it should be

- Every genbank file, that has "." in their name or accession, was changed with underscore to prevent potential problems with input information (Some programs aren't tolerable with dots). 
