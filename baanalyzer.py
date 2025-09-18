from typing import Union, Optional, Dict
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition, Mol
import copy
import json
from .baring import BARing

import pandas as pd
from collections import OrderedDict
import os

consider_ab_site: tuple = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19)


def generate_site_idx_name_str(site_indices: list[int], sub_name: str) -> str:
    """Convert site indices and substituent name to a standardized string format.

    Args:
        site_indices: List of 0-based atom indices
        sub_name: Name of the substituent

    Returns:
        str: Formatted string in "site1,site2,..._subname" format (1-based indices)
    """
    # Convert 0-based indices to 1-based and join with commas
    sites_str = ','.join(str(idx + 1) for idx in site_indices)
    return f"{sites_str}_{sub_name}"


class BAAnalyzer:
    def __init__(self, BA_mol):
        self.BA_mol = BA_mol
        self.BA_type: int = 0
        self.BA_c_atom_num: int = 0
        self.standardized_BA_mol = None
        self.standardized_BA_SMILES: str = ''
        self.Is_standard: bool = False
        self.has_unknown_modification: bool = False
        self.modification_info: dict = {}
        self.cip_tag_info: dict = {}
        self.anomer_tag_info = {}
        self.double_bond_dict = {}
        self.double_bond_atom_list = []
        self.double_bond_number: int = -1
        self._modification_reference_dict: dict = {}
        self.BA_templates: OrderedDict = OrderedDict()
        self._matched_BA_template: dict = {}
        self._CW_map_dict: Dict = {}
        self._init_func()

    def _init_func(self):
        current_dir = os.path.dirname(__file__)
        with open(os.path.join(current_dir, 'conf/FindModification.json'), mode='r') as f:
            self._modification_reference_dict = json.load(f)
        with open(os.path.join(current_dir, 'conf/CW.json'), mode='r') as f:
            self._CW_map_dict = json.load(f)
        with open(os.path.join(current_dir, 'conf/StandardBondIdx.json'), mode='r') as f:
            self._double_bond_idx = json.load(f)
        self._generate_BA_template_mol()

    def analyze(self, BAType: int = None, verbose: bool = False):
        self.standardization_BA_mol(BAType=BAType)
        if self.Is_standard:
            self.find_all_modifications()
            self.find_all_chiral_centers()
            self.find_all_double_bond()
            self.find_all_anomers_tag()
            if verbose:
                self.BADebug()
            return True
        else:
            return False

    def _generate_BA_template_mol(self) -> None:
        for c_atom_length in range(35, 22, -1):  # 有序的 c个数从大到小 类型从1到3
            for BAType in [1, 2, 3]:
                self.BA_templates[(c_atom_length, BAType)] = (
                    BARing.Generate_BAMol_Core(bone_length=c_atom_length, BAType=BAType)
                )

    def _preprocessing_BAMol(self):
        self.BA_mol = Chem.RemoveHs(self.BA_mol)

    def standardization_BA_mol(self, BAType: int = None) -> Mol | None:
        """
        Standardizes atom indexing order by core template matching.

        Args:
            BAType: Specific molecule type to try (1/2/3). If None, tries all types.

        Returns:
            Chem.Mol: Standardized molecule if successful, None if no match found.

        Notes:
            1. Sets StandardIdx property on core atoms (1-based)
            2. Updates isStandard, CAtomNumber, _CoreBAMol attributes
            3. For General type validates terminal atom total degree = 4
        """
        self._preprocessing_BAMol()
        # Validate molecule existence
        if not hasattr(self, 'BA_mol') or self.BA_mol is None:
            raise ValueError("No RDKit molecule found for standardization.")

        # Determine which types to try
        BA_types = [BAType] if BAType is not None else [1, 2, 3]  # 强制使用类型

        # Find all template matches
        type_matches = {}
        for (c_len, b_type), template in self.BA_templates.items():
            if b_type not in BA_types:
                continue
            matches = self.BA_mol.GetSubstructMatches(template)
            if matches:
                type_matches[(c_len, b_type)] = matches[0]

        # Process if matches found
        if not type_matches:
            return None

        try:
            # Get the match with maximum carbon length
            # Assuming we want the template with largest core size
            max_core_key = max(type_matches.keys(), key=lambda x: x[0])  #只比较碳链
            core_atom_indices = type_matches[max_core_key]
            max_c_len, min_BA_type = max_core_key

            # Get side chain indices (excluding core atoms)
            side_chain_indices = [
                atom.GetIdx()
                for atom in self.BA_mol.GetAtoms()
                if atom.GetIdx() not in core_atom_indices
            ]

            # Renumber atoms: core first, then side chains
            new_order = list(core_atom_indices) + side_chain_indices
            new_mol = Chem.RenumberAtoms(self.BA_mol, new_order)

            # Update object attributes
            self.standardized_BA_mol = new_mol
            self.Is_standard = True
            self.BA_c_atom_num = max_c_len
            self.BA_type = min_BA_type
            self._matched_BA_template = self.BA_templates[(max_c_len, min_BA_type)]
            self.standardized_BA_SMILES = Chem.MolToSmiles(new_mol)
            return new_mol

        except (KeyError, IndexError, ValueError) as e:
            Warning(f"Standardization failed: {e}")
            return None

    def get_BA_type(self):
        if self.BA_type == 0:
            return False
        else:
            return self.BA_type

    def find_all_modifications(self, use_chirality: bool = True) -> dict[str, dict]:
        """Identify and record all modifications on the molecule using R-group decomposition,
        including counts for each modification type at each site.

        Args:
            use_chirality: Whether to consider chirality in substructure matching

        Returns:
            A dictionary where:
            - Key is the site index (0-based)
            - Value is another dictionary where:
                - Key is the modification type
                - Value is the count of that modification at the site
        """
        # Perform R-group decomposition
        rgd, fails = rdRGroupDecomposition.RGroupDecompose([self._matched_BA_template], [self.standardized_BA_mol])
        core_mol = next(iter(rgd[0].values()))  # Get the core molecule

        # Step 1: Identify all R-group positions in the core
        core_info = {}
        for atom in core_mol.GetAtoms():
            if 'R' in atom.GetSymbol():
                if len(atom.GetNeighbors()) > 1:
                    raise ValueError("R-group atom has more than one neighbor")

                neighbor = atom.GetNeighbors()[0]
                core_info[atom.GetSymbol()] = {
                    'DummyIdx': atom.GetIdx(),
                    'BASite': neighbor.GetIdx()
                }

        # Step 2: Process each R-group
        modification_dict = {}  # Will store site -> {mod_type: count}
        for R_group_label, R_group_mol in list(rgd[0].items())[1:]:
            # Clean up the R-group molecule
            for atom in R_group_mol.GetAtoms():
                atom.SetAtomMapNum(0)
            R_group_mol = Chem.RemoveHs(R_group_mol)
            Chem.SanitizeMol(R_group_mol)

            # Split into individual substituents if needed
            smiles_parts = Chem.MolToSmiles(R_group_mol).split('.')
            for smiles in smiles_parts:
                # Identify the substituent type
                mol = Chem.MolFromSmiles(smiles)
                match_result = self._R_with_dummy_match_sub(mol, use_chirality)

                # Get the matching template
                try:
                    template = self._modification_reference_dict[match_result['sub']]
                except KeyError:
                    template = match_result['sub']

                # Find matches in the R-group molecule
                pattern = Chem.MolFromSmarts(template)
                matches = R_group_mol.GetSubstructMatches(pattern)

                # Process each match
                site_mod_counts = {}  # Temporary storage for counts per site
                for match in matches:
                    # Find core attachment sites for this match
                    core_sites = []
                    for atom_idx in match:
                        atom = R_group_mol.GetAtomWithIdx(atom_idx)
                        if 'R' in atom.GetSymbol():
                            core_sites.append(core_info[atom.GetSymbol()]['BASite'])
                    if core_sites:
                        site_key = tuple(sorted(core_sites))  # Use tuple of sites as key
                        if site_key not in site_mod_counts:
                            site_mod_counts[site_key] = {}
                        mod_type = match_result['sub']
                        site_mod_counts[site_key][mod_type] = site_mod_counts[site_key].get(mod_type, 0) + 1
                # print(site_mod_counts)
                for all_site, sub_num in site_mod_counts.items():
                    for site in all_site:
                        if site not in self.modification_info:
                            self.modification_info[site] = sub_num
                            self.modification_info[site].update(sub_num)
        return modification_dict

    def find_all_chiral_centers(self) -> Dict[int, str]:
        """Identify and record all chiral centers in the molecule.

        Returns:
            Dictionary mapping atom indices (0-based) to their CIP codes ('R'/'S')

        Notes:
            - Uses RDKit's 0-based atom indexing
            - Clears any existing CIP codes before reassignment
            - Stores results in both cipTagDict and modificationRecorder
        """
        # Clear any existing CIP codes
        for atom in self.standardized_BA_mol.GetAtoms():
            atom.ClearProp("_CIPCode")

        # Find and store all chiral centers
        self.cip_tag_info.clear()  # Reset the dictionary
        for atom_idx, cip_code in Chem.FindMolChiralCenters(self.standardized_BA_mol, includeUnassigned=True):
            self.cip_tag_info[atom_idx] = cip_code
            # self.modificationRecorder.RecordSiteRS(site=atom_idx, RSLabel=cip_code)
        return self.cip_tag_info

    def find_all_double_bond(self):
        """
        Identifies and records all carbon-carbon double bonds in the molecule, including their E/Z stereochemistry.

        Updates:
        - self.doubleBondDict: Dictionary mapping bond atom indices to E/Z labels
        - self.doubleBondAtomList: List of atom indices involved in double bonds (1-based)
        - self.CCDoubleBondNumber: Count of carbon-carbon double bonds
        - modificationRecorder: Records each double bond with its stereochemistry
        """

        for bond in self.standardized_BA_mol.GetBonds():
            # Check if this is a C=C double bond
            if (bond.GetBondTypeAsDouble() > 1 and
                    bond.GetBeginAtom().GetSymbol() == 'C' and
                    bond.GetEndAtom().GetSymbol() == 'C'):

                idx1, idx2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                idx1, idx2 = min(idx1, idx2), max(idx1, idx2)
                bond_key = (idx1, idx2)
                # Initialize with None stereochemistry
                self.double_bond_dict[bond_key] = None
                self.double_bond_atom_list.extend([idx1 + 1, idx2 + 1])  # Convert to 1-based
                if self.double_bond_atom_list is not None:
                    # Determine stereochemistry
                    ez_label = None
                    if bond.GetStereo() == Chem.BondStereo.STEREOE:
                        ez_label = 'E'
                        self.double_bond_dict[bond_key] = 'E'
                    elif bond.GetStereo() == Chem.BondStereo.STEREOZ:
                        ez_label = 'Z'
                        self.double_bond_dict[bond_key] = 'Z'
        self.double_bond_number = len(self.double_bond_dict)

    def find_all_anomers_tag(self):
        for site in self.cip_tag_info.keys():
            if site not in consider_ab_site:
                continue
            if self.cip_tag_info[site] == '?':
                self.anomer_tag_info[site] = '?'
                continue
            self.anomer_tag_info[site] = self.cal_site_anomers_tag(site)

    def cal_site_anomers_tag(self, site: int, mol: object = None) -> Optional[str]:
        """Determine the anomeric configuration (alpha/beta) for a given site in a bile acid molecule.

        Args:
            site: 0-based atom index to analyze
            mol: Optional molecule to analyze (uses self.BAMol if None)

        Returns:
            'a' for alpha, 'b' for beta, or None if not applicable

        Notes:
            - Works with standardized bile acid molecules
            - Uses CIP rules to determine stereochemistry
            - Handles two modes: sites with H (mode 1) and without H (mode 2)
        """
        # Ensure chiral centers are identified
        self.find_all_chiral_centers()

        # Prepare the molecule
        working_mol = copy.deepcopy(mol if mol else self.standardized_BA_mol)
        working_mol = Chem.AddHs(working_mol)

        # Define site classifications
        ANOMER_SITES_WITH_H = {0, 1, 2, 3, 5, 6, 10, 11, 14, 15, 16, 17, 18, 20, 21, 22}
        ANOMER_SITES_WITHOUT_H = {4, 7, 8, 9, 12, 13, 16, 19}  # lipidmaps substitute sites

        # Validate input site
        if site not in ANOMER_SITES_WITH_H | ANOMER_SITES_WITHOUT_H:
            return None
        if site not in self.cip_tag_info:
            return None

        # Determine analysis mode
        mode = 1 if site in ANOMER_SITES_WITH_H else 2

        # Prepare standard BA atom list (excluding special cases)
        standard_ba_atoms = set(range(self.BA_c_atom_num)) - {17, 18}

        # Assign CIP ranks
        Chem.AssignStereochemistry(working_mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

        # Analyze neighbors
        center_atom = working_mol.GetAtomWithIdx(site)
        neighbor_ranks = {}

        for neighbor in center_atom.GetNeighbors():
            idx = neighbor.GetIdx()
            cip_rank = int(neighbor.GetProp('_CIPRank'))
            key = 'x' if idx not in standard_ba_atoms else str(idx)
            neighbor_ranks[key] = cip_rank

        # Sort neighbors by CIP rank (descending)
        sorted_neighbors = sorted(neighbor_ranks.items(), key=lambda x: x[1], reverse=True)
        neighbor_str = '_'.join(k for k, v in sorted_neighbors)

        # Determine configuration
        if mode == 1:
            # Mode 1: Site has hydrogen (H must be lowest CIP rank)
            is_clockwise = neighbor_str in self._CW_map_dict[str(site)]
            neighbor_wise = 'S' if is_clockwise else 'N'
            rs_rule = 'S' if self.cip_tag_info[site] == 'R' else 'N'
            return 'b' if rs_rule == neighbor_wise else 'a'

        else:  # mode == 2
            # Mode 2: Site without hydrogen
            *main_neighbors, min_neighbor = neighbor_str.split('_')
            main_neighbor_str = '_'.join(main_neighbors)

            # Get possible clockwise arrangements for this minimum rank
            cw_options = self._CW_map_dict[str(site)].get(min_neighbor, [])
            is_clockwise = main_neighbor_str in cw_options
            neighbor_wise = 'S' if is_clockwise else 'N'
            rs_rule = 'S' if self.cip_tag_info[site] == 'R' else 'N'

            return 'b' if rs_rule == neighbor_wise else 'a'

    def get_standardized_BA_mol(self):
        return self.standardized_BA_mol

    def _R_with_dummy_match_sub(self, mol, use_chirality) -> dict:
        """Find matching modification templates and count dummy atoms in the molecule.

        Args:
            mol: RDKit molecule object
            use_chirality: bool, whether to use chirality in substructure matching

        Returns:
            dict: {'dummyNumber': count of dummy atoms, 'sub': matched substructure SMILES}
            :param mol:
            :param use_chirality:
            :return:
        """
        match_result = {'dummyNumber': 0, 'sub': ''}

        # First pass: try to match known templates
        for smarts in self._modification_reference_dict.values():
            template = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(template, useChirality=use_chirality)
            # We want exact single match (full molecule match)
            if len(matches) == 1 and len(matches[0]) == mol.GetNumAtoms():
                dummy_count = sum(1 for atom_idx in matches[0]
                                  if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 0)

                if dummy_count > 0:  # Only store if we found dummies
                    match_result.update({
                        'dummyNumber': dummy_count,
                        'sub': smarts
                    })
                    return match_result  # Return first valid match

            # Second pass: count dummies in the whole molecule if no template matched
            dummy_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0)
            if dummy_count > 0:
                match_result.update({
                    'dummyNumber': dummy_count,
                    # 'sub': Chem.MolToSmarts(mol)
                    'sub': smarts
                })
                self.has_unknown_modification = True
        return match_result

    def update_mol(self, new_mol):
        self.BA_mol = new_mol
        return True

    def debug(self):
        if not hasattr(self, 'modification_info') or not self.modification_info:
            print("No modification information available")
            return 0
        print(f'BA mol SMILES: {Chem.MolToSmiles(self.standardized_BA_mol)}')
        print(f'BA Type: {self.BA_type}')
        for sub, site_num in self.modification_info.items():
            for site_num_ in site_num:
                site = site_num_[0]
                num = int(site_num_[1])
                print(f'Sub : {sub},site: {site},num:{num}')
        print(f'DB Num {self.double_bond_number}')
        print(f'DB Info {self.double_bond_dict}')
        print(f'AB Info {self.anomer_tag_info}')
        print(f'C Num {self.BA_c_atom_num}')
        return 0

    def get_site_RS_tag(self, site: int, start_index: int = 1) -> str:
        """

        :param start_index: 0 or 1
        :param site: site number
        :return: str R/S
        """
        if start_index:
            site = site - 1
        try:
            RS_tag = self.cip_tag_info[site]
        except KeyError:
            RS_tag = 'None'
        return RS_tag

    def get_site_AB_tag(self, site: int, start_index: int = 1) -> str:
        """

        :param start_index:
        :param site: site number from 0
        :return: str R/S
        """
        if start_index:
            site = site - 1
        AB_tag = ''
        try:
            AB_tag = self.anomer_tag_info[site]
        except KeyError:
            AB_tag = 'None'
        return AB_tag

    def change_mod_name_in_abbr_ref(self):
        pass

    def get_site_mod(self, site: int, start_index: int = 1) -> str:
        """

        :param start_index:
        :param site: site number from 0
        :return: str R/S
        """
        if start_index:
            site = site - 1
        try:
            mod_str = self.modification_info[site]
        except KeyError:
            return 'H'
        swapped_dict = {value: key for key, value in self._modification_reference_dict.items()}
        mod_str_list = []
        for mod in mod_str:
            try:
                mod_name = swapped_dict[mod]
            except KeyError:
                # print(f"Warning: Unknown modification reference '{mod_str}'")
                mod_name = mod_str
            if mod_str[mod] > 1:
                mod_name = f'{mod_str[mod]}-{mod_name}'
            else:
                mod_name = mod_name
            mod_str_list.append(mod_name)
        return ','.join(mod_str_list)

    def BADebug(self) -> bool:
        """
        Debugging utility to print modification information for the molecular structure.

        Prints each modification site with its stereo chemical configuration (RS/AB tags)
        and modification name in a formatted way.

        Returns:
            bool: False if the structure is not standard, True otherwise.
        """
        print(f'SMILES:{self.standardized_BA_SMILES}')
        print(f'BA Type: {self.BA_type}')
        if self.double_bond_dict:
            print(f'Double bond: {list(self.double_bond_dict.keys())}')
        for site in range(0, self.BA_c_atom_num + 1):
            self.get_site_mod(site, start_index=0)
            print(f'Site:{site + 1:<5} RS:{self.get_site_RS_tag(site, start_index=0):<5} '
                  f'AB:{self.get_site_AB_tag(site, start_index=0):<5} Mod: {self.get_site_mod(site, start_index=0):<5}')
        return True

    def get_BA_mol_info_dataframe(self, stero_tag_type: str = 'AB') -> pd.DataFrame:
        """
        Generate a molecular information dataframe for BA (Bile Acid) molecules.

        Args:
            stero_tag_type (str): Stereo labeling type, either 'AB' or 'RS'. Defaults to 'AB'.

        Returns:
            pd.DataFrame: A dataframe containing molecular information with stereo tags and modifications.

        Raises:
            ValueError: If an invalid stereo_tag_type is provided.
        """
        # Initialize dataframe columns
        columns = ['Double Bond', 'BA Type']
        carbon_atoms = [str(i) for i in range(1, self.BA_c_atom_num + 1)]
        columns.extend(carbon_atoms)

        # Create empty dataframe with specified columns
        info_dataframe = pd.DataFrame(columns=columns, index=['SMILES'])

        # Validate stereo tag type
        if stero_tag_type not in ['AB', 'RS']:
            raise ValueError('Stereo_tag_type must be either "AB" or "RS"')

        # Determine which sites to consider based on tag type
        considered_sites = consider_ab_site if stero_tag_type == 'AB' else range(1, self.BA_c_atom_num + 1)

        # Populate dataframe with site information
        for site in range(1, self.BA_c_atom_num + 1):
            if site in considered_sites:
                # Get appropriate stereo tag
                if stero_tag_type == 'AB':
                    tag = self.get_site_AB_tag(site, start_index=1)
                else:
                    tag = self.get_site_RS_tag(site, start_index=1)

                # Format the site information
                site_mod = self.get_site_mod(site, start_index=1)
                info_dataframe.loc['SMILES', str(site)] = f'{tag}-{site_mod}'

        return info_dataframe


def _change_H_to_dummy(mol, site, order) -> Union[bool, Chem.Mol]:
    """Replace hydrogen at specified site with a dummy atom with given bond order.

    Args:
        mol: RDKit molecule object
        site: int, index of atom to modify
        order: int (1-3), desired bond order to dummy atom

    Returns:
        RWMol: Modified molecule or False if valence exceeds
    """
    # Validate input
    if order not in {1, 2, 3}:
        raise ValueError("Order must be 1, 2, or 3")

    # Create working copy
    mol_copy = Chem.RWMol(mol)
    periodic_table = Chem.GetPeriodicTable()
    site_atom = mol_copy.GetAtomWithIdx(site)

    # Check valence
    default_valence = periodic_table.GetDefaultValence(site_atom.GetSymbol())
    if site_atom.GetDegree() >= default_valence:
        print('Valence will exceed default')
        return False

    # Prepare the site by adding H if needed
    mol_copy = Chem.AddHs(Chem.RWMol(mol_copy), onlyOnAtoms=[site])

    # Find the hydrogen to replace
    h_neighbor = next((a for a in site_atom.GetNeighbors()
                       if a.GetSymbol() == 'H'), None)
    if not h_neighbor:
        print('No hydrogen found at site')
        return False

    h_idx = h_neighbor.GetIdx()

    # Replace H with dummy atom
    dummy = Chem.Atom(0)
    dummy.SetAtomMapNum(1)
    mol_copy.ReplaceAtom(h_idx, dummy)

    # Adjust bond order and hydrogens based on requested order
    if order == 1:
        pass  # Default is single bond
    else:
        bond = mol_copy.GetBondBetweenAtoms(h_idx, site)
        bond.SetBondType(Chem.BondType.DOUBLE if order == 2 else Chem.BondType.TRIPLE)
        site_atom.SetNumExplicitHs(0)
        site_atom.SetNoImplicit(False)
        mol_copy = Chem.RemoveHs(mol_copy)
    return mol_copy
