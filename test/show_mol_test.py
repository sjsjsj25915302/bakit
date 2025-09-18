from utils import plotMol
from rdkit import Chem

SMILES = 'C[C@@H]([C@H]1CC[C@@H]2[C@]1(C)CC[C@H]3[C@H]2[C@H](O)C[C@H]4[C@]3(C)CC[C@@H](O)C4)CCC(O[C@H]5[C@H](O)[' \
         'C@@H](O)[C@H](O)[C@@H](C(O)=O)O5)=O'

mol = Chem.MolFromSmiles(SMILES)

mol_draw = plotMol.MolDraw(mol)
mol_draw.mol_plot()

mol_draw = plotMol.MolDraw(mol, rotate=98)
mol_draw.mol_plot()

mol_draw = plotMol.MolDraw(mol, rotate=98, show_atom_index=True,
                           show_bond_index=True, show_stereo_label=True)
mol_draw.mol_plot()


mol_draw = plotMol.MolDraw(mol, rotate=98, save_path='fig/plot_mol_test_1')
mol_draw.mol_save()

mol_draw = plotMol.MolDraw(mol, rotate=98, high_light_atoms=[1,2,3,4,5,6],
                           high_light_bonds=[1,2,3,4,5,6]
                           )

mol_draw.mol_plot()

mol_draw = plotMol.MolDraw(mol, size=(200,200), use_acs=True, rotate=98, show_atom_index=True
                           , save_path='fig/plot_mol_test_acs')
mol_draw.mol_plot()
mol_draw.mol_save()