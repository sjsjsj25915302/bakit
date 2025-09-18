"""
构建BA分子环模板
# encoding: utf-8
# by SunJian
# version 0.0.2
"""
from rdkit.Chem.rdchem import RWMol
from rdkit import Chem


class BAMolType:
    General: int = 1
    NoOH: int = 2
    NoOOH: int = 3


class BARing:
    def __init__(self):
        pass

    @staticmethod
    def Generate_BAMol_Core(bone_length: int, BAType: int = False) -> object:
        """
            First, a ring skeleton's chain form is generated, and then atoms are sequentially add onto the skeleton.
            :return mol object
        """
        assert bone_length >= 22, 'Length must longer than 22'
        BA_Mol = RWMol()
        for _ in range(0, 17):
            BA_Mol.AddAtom(Chem.Atom(6))
        for i in range(0, 16):
            BA_Mol.AddBond(i, i + 1, Chem.BondType.UNSPECIFIED)
        remaining_carbons = bone_length - 17
        for _ in range(0, remaining_carbons):
            BA_Mol.AddAtom(Chem.Atom(6))
        ring_bonds = [(0, 9), (4, 9), (7, 13), (8, 10), (12, 16), (16, 19), (19, 20),
                      (19, 21), (9, 18), (12, 17)]
        for bond in ring_bonds:  # connect ring
            BA_Mol.AddBond(*bond, Chem.BondType.SINGLE)
        BA_Mol.RemoveBond(9, 10)
        for i in range(21, bone_length - 1):  # connect chain
            BA_Mol.AddBond(i, i + 1, Chem.BondType.UNSPECIFIED)

        for bond in BA_Mol.GetBonds():
            bond.SetBondType(Chem.BondType.UNSPECIFIED)


        if BAType == BAMolType.General:
            newIdx = BA_Mol.AddAtom(Chem.Atom(8))
            BA_Mol.AddBond(beginAtomIdx=newIdx, endAtomIdx=bone_length - 1, order=Chem.BondType.DOUBLE)
            newIdx = BA_Mol.AddAtom(Chem.Atom(8))
            BA_Mol.AddBond(beginAtomIdx=newIdx, endAtomIdx=bone_length - 1, order=Chem.BondType.SINGLE)
            BA_Mol.UpdatePropertyCache()

        if BAType == BAMolType.NoOH:
            newIdx = BA_Mol.AddAtom(Chem.Atom(8))
            BA_Mol.AddBond(beginAtomIdx=newIdx, endAtomIdx=bone_length - 1, order=Chem.BondType.DOUBLE)
            BA_Mol.UpdatePropertyCache()

        if BAType == BAMolType.NoOOH:
            BA_Mol.UpdatePropertyCache()

        Chem.SanitizeMol(BA_Mol)
        return BA_Mol.GetMol()
