"""
test file for BAAnalyzer standardization and analyzer
"""
from bakit import BAAnalyzer, MolDraw
from rdkit import Chem
import pandas as pd

BA_SMILES_list_pd = pd.read_csv('./test/test.csv')
BA_SMILES_list = BA_SMILES_list_pd['SMILES'].to_list()
BA_SMILES_list = ['C[C@@H]([C@H]1CC[C@H]2C3=CC[C@H]4C[C@H](CC[C@@]4([C@H]3CC[C@]12C)C)O)CCC[C@@H](C(O)=O)C']
for SMILES in BA_SMILES_list:
    test_mol = Chem.MolFromSmiles(SMILES)
    BAMol = BAAnalyzer(test_mol)
    BAMol.get_site_mod(24,start_index=1)
    if BAMol.analyze(verbose=True):
        mol_draw = MolDraw(show_atom_index=True, show_stereo_label=True, cairo=True,
                           title='ST', rotate=137)

        print(BAMol.get_BA_mol_info_dataframe())
        mol_draw.mol_draw(BAMol.get_standardized_BA_mol(), )
    else:

        print(f'Failed {SMILES}')
