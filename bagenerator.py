from .baanalyzer import BAAnalyzer
from rdkit import Chem
from rdkit.Chem.rdchem import RWMol
import copy
from typing import Optional
import json


class BAGenerator(BAAnalyzer):
    def __init__(self, BA_mol: Optional[Chem.Mol] = None, BAType:  Optional[int] = None) -> None:
        """Initialize BAGenerator with an optional BA molecule."""
        super().__init__(BA_mol)
        self.initial_BA_mol = copy.deepcopy(self.BA_mol)
        self.standardization_BA_mol(BAType=BAType)
        self.analyze(BAType=BAType,verbose=False)
        self.final_BA_mol: Chem.RWMol = copy.deepcopy(self.standardized_BA_mol)

        with open('src/bakit/conf/AddModification.json', mode='r') as f:
            self._add_mod_ref_dict = json.load(f)

    def denovo_generate_BA_mol(self, BA_smiles: str) -> str:
        """Placeholder method for de novo BA molecule generation."""
        pass

    def get_core_BA_mol_from_SMILES(self, BA_smiles: str) -> bool:
        """
        Generate core BA molecule from SMILES string and standardize it.

        Args:
            BA_smiles: SMILES string representation of the BA molecule

        Returns:
            bool: True if standardization successful, False otherwise
        """
        self.update_mol(Chem.MolFromSmiles(BA_smiles))
        self.standardization_BA_mol()
        return self.Is_standard

    def get_core_ba_mol_from_mol(self, BA_mol: Chem.Mol) -> bool:
        """
        Generate core BA molecule from RDKit Mol object and standardize it.

        Args:
            BA_mol: RDKit Mol object representing the BA molecule

        Returns:
            bool: True if standardization successful, False otherwise
        """
        if not isinstance(BA_mol, Chem.Mol):
            return False
        self.update_mol(BA_mol)
        self.standardization_BA_mol()
        return self.Is_standard

    def change_chiral_centers(self, site: int, start_index: int = 0) -> bool:
        """
        Invert chirality at a specific chiral center site.

        Args:
            site: Atom index of the chiral center
            start_index: Index offset for site numbering (default: 0)

        Returns:
            bool: True if inversion successful, False otherwise
        """
        adjusted_site = site - start_index
        # BAMol_ = copy.deepcopy(self.BA_mol)
        atom = self.BA_mol.GetAtomWithIdx(adjusted_site)

        # Clear CIP code to allow recalculation
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")

        # Invert chirality if atom is chiral
        if atom.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                                   Chem.ChiralType.CHI_TETRAHEDRAL_CW):
            atom.InvertChirality()
            self.update_mol(self.BA_mol)
            self.analyze()
            return True

        return False

    def change_site_RS_tag(self, site: int, RS_tag: str, start_index: int = 0) -> bool:
        """
        Change the R/S tag for a specific chiral center site.

        Args:
            site: Atom index of the chiral center
            RS_tag: Desired chirality tag ('R' or 'S')
            start_index: Index offset for site numbering (default: 0)

        Returns:
            bool: True if operation successful, False otherwise
        """
        adjusted_site = site - start_index
        # BAMol_ = copy.deepcopy(self.BA_mol)
        atom = self.standardized_BA_mol.GetAtomWithIdx(adjusted_site)

        # Validate inputs
        if adjusted_site not in self.cip_tag_info:
            print('No CIP information for this site')
            return False

        if RS_tag not in ['R', 'S']:
            return False

        # Clear existing CIP code to allow recalculation
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")

        # Set initial chirality and update molecule
        atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)


        # Check if desired RS tag was achieved
        if RS_tag == self.get_site_RS_tag(adjusted_site):
            pass
        else:
            # Invert chirality if needed
            self.change_chiral_centers(adjusted_site)
        try:
            self.update_mol(self.standardized_BA_mol)
            self.analyze()
            self.final_BA_mol = self.standardized_BA_mol
        except Exception as e:
            print(f'Add Modification Failed: {str(e)}')
            return False

    def change_site_AB_tag(self, site: int, ab_tag: str, start_index: int = 0) -> bool:
        """
        Change the anomeric (a/b) tag for a specific site.

        Args:
            site: Atom index of the anomeric center
            ab_tag: Desired anomeric configuration ('a' or 'b')
            start_index: Index offset for site numbering (default: 0)

        Returns:
            bool: True if operation successful, False otherwise
        """
        adjusted_site = site - start_index
        # BAMol_ = copy.deepcopy(self.BA_mol)
        atom = self.standardized_BA_mol.GetAtomWithIdx(adjusted_site)

        # Validate anomeric site and tag
        if adjusted_site not in self.anomer_tag_info:
            return False

        if ab_tag not in ['a', 'b']:
            return False

        # Clear CIP code to allow recalculation if present
        if atom.HasProp("_CIPCode"):
            atom.ClearProp("_CIPCode")

        # Handle unspecified chirality
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        try:
            self.update_mol(self.standardized_BA_mol)
            self.analyze()
            self.final_BA_mol = self.standardized_BA_mol
        except Exception as e:
            print(f'Add Modification Failed: {str(e)}')
            return False

        # Check if desired configuration is already present
        if ab_tag == self.get_site_AB_tag(adjusted_site):
            return True
        else:
            # Invert chirality to achieve desired configuration
            self.change_chiral_centers(adjusted_site)
            return True

    def change_site_mod(self, site: int, mod: str, start_index: int = 0) -> bool:
        """
        Modify specific sites on the molecule with the given modification.

        Args:
            site: Site index or tuple of indices to modify (1-based)
            mod: Abbreviation of the modification to apply
            start_index: Starting index for modification (default 0)

        Returns:
            bool: True if modification was successful, False otherwise
        """
        # BAMol_ = copy.deepcopy(self.BA_mol)
        # mol_draw = plotMol.MolDraw(self.standardized_BA_mol, show_atom_index=True, show_stereo_label=True, title='1',
        #                            rotate=137)
        # mol_draw.mol_plot()
        # Get modification reference pattern
        try:
            add_mod_smart = self._add_mod_ref_dict[mod]
        except KeyError:
            print(f'No modification reference for {mod}')
            return False
        add_mod_mol = Chem.MolFromSmiles(add_mod_smart)

        # Process site indices (convert to 0-based)
        if isinstance(site, int):
            site_list = [site - start_index]
        elif isinstance(site, tuple):
            site_list = [x - start_index for x in site]
        else:
            print(f'Invalid site type: {type(site)}')
            return False

        # Add hydrogens only to the target atoms
        self.final_BA_mol = Chem.AddHs(self.final_BA_mol, onlyOnAtoms=site_list)
        BA_RWMol = RWMol(self.final_BA_mol)
        # Process each modification site
        for i, site_idx in enumerate(site_list):
            # Find valence from substituent molecule
            valence = None
            for atom in add_mod_mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    valence = atom.GetTotalValence()
                    break

            if valence in (1, 2):
                # Find and replace hydrogen atom with dummy atom
                toDummyHAtomSite = None
                for atom in BA_RWMol.GetAtomWithIdx(site_idx).GetNeighbors():
                    if atom.GetSymbol() == 'H':
                        toDummyHAtomSite = atom.GetIdx()
                        break

                if toDummyHAtomSite is None:
                    print(f'No hydrogen atom found at site {site_idx + 1}')
                    return False

                BA_RWMol.ReplaceAtom(toDummyHAtomSite, Chem.Atom(0))
                mol_ = BA_RWMol.GetMol()
                mol_.GetAtomWithIdx(toDummyHAtomSite).SetAtomMapNum(i + 1 + start_index)
                BA_RWMol = RWMol(mol_)

                # Adjust hydrogen count and implicit properties
                BA_RWMol.GetAtomWithIdx(site_idx).SetNumExplicitHs(0)
                BA_RWMol.GetAtomWithIdx(site_idx).SetNoImplicit(False)
                BA_RWMol = Chem.RemoveHs(BA_RWMol)

        # Combine molecules and apply modification
        combineMol = Chem.CombineMols(add_mod_mol, BA_RWMol)
        if combineMol:
            combineMol = Chem.RemoveHs(combineMol)
            combineMol = Chem.molzip(combineMol)
        # mol_draw = plotMol.MolDraw(combineMol, show_atom_index=True, show_stereo_label=True,
        #                            rotate=137)
        # mol_draw.mol_plot()
        # Verify the modification was successful
        try:
            self.update_mol(combineMol)
            self.analyze()
            self.final_BA_mol = self.standardized_BA_mol
        except Exception as e:
            print(f'Add Modification Failed: {str(e)}')
            return False

        return True


    def get_final_BA_mol(self):
        return self.standardized_BA_mol

    def change_BA_type(self,BA_type :int) -> bool:
        pass

