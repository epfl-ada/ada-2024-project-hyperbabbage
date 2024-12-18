from os.path import join
DATA_RAW = 'data/raw/'
DATA_CLEAN = 'data/clean/'

BINDINGDB_RAW = DATA_RAW + 'BindingDB_All.tsv'
BINDINGDB_CLEAN = DATA_CLEAN + 'BindingDB_Cleaned.pkl'

DRUGBANK_XML = DATA_RAW + 'full_database.xml'
DRUGBANK_LIGAND_PARSED = DATA_CLEAN + 'parsed_DrugBank_ligand.pkl'
DRUGBANK_PROTEIN_PARSED = DATA_CLEAN + 'parsed_DrugBank_protein.pkl'
LIGANDS_RELATED_TO_PROTEIN = DATA_CLEAN + 'ligands_related_to_cancer_proteins.pkl'

MERGED = DATA_CLEAN + 'merged_dataframe.pkl'
COSMIC_PROTEINS = DATA_RAW + 'protein_class_COSMIC.tsv'

DOI_DF_PATH = DATA_CLEAN +'df_doi.pkl'
