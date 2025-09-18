"""
test file for BARing
"""
from bakit import MolDraw
from bakit import BARing
for i in [1,2,3]:
    draw = MolDraw(rotate=180,cairo=True)
    draw.mol_draw(mol=BARing.Generate_BAMol_Core(bone_length=24, BAType=i))




