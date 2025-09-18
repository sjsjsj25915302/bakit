"""
test file for BAGenerator
"""
from bakit import BAGenerator
from rdkit import Chem
from bakit import plotMol

core_mol = 'C[C@H](CCC(O)=O)[C@H]1CC[C@H]2[C@@H]3CC[C@@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C'
BA_mol = Chem.MolFromSmiles(core_mol)

mol_draw = plotMol.MolDraw(show_atom_index=True, show_stereo_label=True,
                           rotate=137, title='core mol')
mol_draw.mol_draw(BA_mol, )
generator = BAGenerator(BA_mol)
mol_draw = plotMol.MolDraw(show_atom_index=True, show_stereo_label=True,
                           rotate=137, title='Before')
mol_draw.mol_draw(generator.BA_mol, )
generator.change_site_mod(site=0, mod='OH')
mol_draw = plotMol.MolDraw(show_atom_index=True, show_stereo_label=True,
                           rotate=137, title='After_mod')
mol_draw.mol_draw(generator.get_standardized_BA_mol(), )
generator.change_site_RS_tag(site=0, RS_tag='R')
mol_draw = plotMol.MolDraw(show_atom_index=True, show_stereo_label=True,
                           rotate=137, title='After_RS')
mol_draw.mol_draw(generator.get_standardized_BA_mol(), )

generator.change_site_AB_tag(site=4, ab_tag='b')
mol_draw = plotMol.MolDraw(show_atom_index=True, show_stereo_label=True,
                           rotate=137, title='After_ab')
mol_draw.mol_draw(generator.get_standardized_BA_mol(), )
