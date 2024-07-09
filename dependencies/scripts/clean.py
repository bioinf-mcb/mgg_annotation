""" remove unnecessaru files """
from pathlib import Path
def cleaning(config):
    
 
    # paths & params
    CLEAN = config['CLEAN']

    OUTPUT_DIR = config['OUTPUT_DIR']
    ORF_PREDICTION_DIR = Path(OUTPUT_DIR, '2_ORF_PREDICTION')
    ANNOTATION_DIR = Path(OUTPUT_DIR, '3_ANNOTATION')
    GENBANK_DIR = Path(OUTPUT_DIR, '4_GENBANK')
    ORFS_DIR = Path(ORF_PREDICTION_DIR, '3_ORFS')
    PROTEINS_DIR = Path(ORF_PREDICTION_DIR, '4_PROTEINS')

    files_to_keep = [Path(GENBANK_DIR),
                     Path(OUTPUT_DIR, 'annotation.tsv'),
                     Path(OUTPUT_DIR, "search.tsv"),
                     Path(OUTPUT_DIR, "report.tsv"),
                     Path(ORFS_DIR),
                     Path(PROTEINS_DIR),
                     Path(ANNOTATION_DIR, "PCs2proteins.tsv"),
                     Path(ANNOTATION_DIR, 'proteins.fasta'),
                     Path(ANNOTATION_DIR, 'msa.a3m'),
                     Path(OUTPUT_DIR, 'ref_info.tsv')]

    # switch
    if CLEAN:
            # run
        print(f'\nCleaning... ', end='')
        delete_files_recursive(Path(OUTPUT_DIR), files_to_keep)
        print('Done!')
    else: pass

    


def delete_files_recursive(folder_path, files_to_keep):
    for item in folder_path.iterdir():
        if item not in files_to_keep:
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                delete_files_recursive(item, files_to_keep)
                try:
                    item.rmdir()
                except OSError as e:
                    # Ignorujemy błąd, jeśli katalog jest nadal niepusty
                    pass

