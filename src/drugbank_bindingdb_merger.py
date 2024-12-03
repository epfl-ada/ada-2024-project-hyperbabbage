from typing import List
import pandas as pd
import os
import numpy as np
from tqdm.notebook import tqdm

"""
The main idea:  We want to maximize the amount of data remaining after merge between BindingDB and DrugBank
                We don't have an identifier that matches them clearly, so we defined a set of identifiers
                If any of those matches between the two databases then we consider them as a match and we merge

Since an immidiate left-join on every identifier will make the computation very slow we merge in two steps:

1)  We execute an inner join on every identifier and append the values under each other
    This will result in some duplicated rows (because a row is highly to match in more than 1 identifier)
2)  We execute a left join between the bindingDB and the merged dataframe => We don't loose data from bindingDB

    a) We need to remove the duplicated rows from point 1)
    b) After left joining we will have duplicated features as well, that should be removed. The duplicated features caused by
    the left-join will be marked with  *_duplicated*, so we know what features to remove
"""
class DrugBank_BindingDB_Merger:

    def __init__(self):
        self.drugbank_df = None
        self.binding_df = None
        self.merged_df = None

    def merge(self, drugbank_df : pd.DataFrame, binding_df : pd.DataFrame):

        self.drugbank_df = drugbank_df
        self.binding_df = binding_df

        #1)
        temp_file = 'temp.csv'
        if os.path.exists(temp_file):
            os.remove(temp_file)

        identifiers = self._rename_cols_and_get_identifiers()
        before_left_merge = self._merge_dataframes_on_identifiers(identifiers, temp_file)
        #2)
        self._left_join(before_left_merge)
        
        return self.merged_df
    
    # Rename columns in BindingDB and in DrugBank to unify naming conventions and return identifiers
    def _rename_cols_and_get_identifiers(self) -> List[str]:
        
        # self.binding_df.rename(columns={
        #     'PubChem CID': 'PubChem_CID',
        #     'ChEBI ID of Ligand': 'ChEBI_ID',
        #     'ChEMBL ID of Ligand': 'ChEMBL_ID',
        #     'DrugBank ID of Ligand': 'DrugBank_ID',
        #     'KEGG ID of Ligand': 'KEGG_ID',
        #     'ZINC ID of Ligand': 'ZINC_ID',
        #     'Ligand SMILES': 'SMILES',
        #     'Ligand InChI Key': 'InChI_Key',
        #     'BindingDB MonomerID': 'BindingDB_ID',
        # }, inplace=True)

        self.drugbank_df.rename(columns={
            'chebi': 'chebi_id',
            'chembl': 'chembl_id',
            'pubchem': 'pubchem_cid',
            'PubChem Substance': 'pubchem_sid',
            'id': 'drugbank_id',
            'bindingdb': 'bindingdb_id',
            'ZINC': 'zinc_id',
            'SMILES': 'smiles',
            'InChI': 'inchi_key',
            'KEGG Compound': 'kegg_id'
        }, inplace=True)


        self.binding_df['Unique_ID'] = np.arange(len(self.binding_df))

        # List of identifiers to merge on
        identifiers = [
            'pubchem_cid',
            'chebi_id',
            'chembl_id',
            'drugbank_id',
            'bindingdb_id',
            # 'zinc_id',
            'smiles',
            # 'inchi',
            'inchi_key'
        ]

        for id in identifiers:
            assert id in self.binding_df.columns, f'{id} not in BindingDB'
            assert id in self.drugbank_df.columns, f'{id} not in DrugBank'
        
        return identifiers

    # Function to process and merge on each identifier individually
    def _merge_dataframes_on_identifiers(self,identifiers, output_file):
        for identifier in tqdm(identifiers):
            if identifier in self.binding_df.columns and identifier in self.drugbank_df.columns:
                
                # Drop rows with NaN in the identifier columns
                binding_df_id = self.binding_df.dropna(subset=[identifier])
                drugbank_df_id = self.drugbank_df.dropna(subset=[identifier])
                
                # Convert identifier columns to string to avoid type mismatches
                binding_df_id[identifier] = binding_df_id[identifier].astype(str)
                drugbank_df_id[identifier] = drugbank_df_id[identifier].astype(str)
                
                # Perform the merge with inner join
                merged_df = pd.merge(
                    binding_df_id, drugbank_df_id, 
                    on=identifier, 
                    how='inner', 
                    suffixes=('_BindingDB', '_DrugBank')
                )
                
                if not merged_df.empty:
                    # Add a column to indicate which identifier was matched
                    merged_df['Matched_On'] = identifier
                    
                    print('run')
                    # Write to CSV in append mode
                    merged_df.to_csv(
                        output_file, 
                        mode='a', 
                        index=False, 
                        header=not os.path.exists(output_file)
                    )
        return_df =  pd.read_csv(output_file)
        return_df.drop_duplicates(subset=['Unique_ID'], inplace=True)
        os.remove(output_file)
        return return_df


    def _left_join(self, merged_df: pd.DataFrame):
        # Left join binding db with merged_df (on Unique_ID), don't duplicate columns though
        binding_readded = pd.merge(self.binding_df, merged_df, on='Unique_ID', how='left', suffixes=('', '_duplicated'))
        all_cols = list(binding_readded.columns)
        cols_to_keep = []

        for col in all_cols:
            if not ('_duplicated' in col and col.split('_duplicated')[0] in all_cols):
                cols_to_keep.append(col)

        # Keep only the columns in cols_to_keep
        self.merged_df = binding_readded[cols_to_keep]