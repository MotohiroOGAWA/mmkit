from __future__ import annotations

from typing import Union, Tuple, Dict
from rdkit import Chem
from rdkit.Chem import ResonanceMolSupplier, ResonanceFlags
from .Formula import Formula
from typing import Optional
from collections import defaultdict
from bidict import bidict

class Compound:
    """
    Represent a chemical compound as an RDKit molecule.

    A ``Compound`` stores a canonicalized RDKit molecule together with
    derived properties such as molecular formula, exact mass, and charge.
    Atom map numbers can also be assigned and tracked.
    """

    def __init__(self, mol: Chem.Mol, overwrite_atom_map: bool = False) -> None:
        """
        Initialize a compound from an RDKit molecule.

        Parameters
        ----------
        mol : Chem.Mol
            RDKit molecule object.
        overwrite_atom_map : bool, default=False
            Whether to overwrite existing atom map numbers.

        Raises
        ------
        AssertionError
            If ``mol`` is not an RDKit molecule.
        ValueError
            If ``mol`` cannot be converted into a valid canonical molecule.
        """
        assert isinstance(mol, Chem.Mol), "mol must be an RDKit Mol object"

        smiles = Chem.MolToSmiles(mol, canonical=True)
        canonical_mol = Chem.MolFromSmiles(smiles)

        if canonical_mol is None:
            raise ValueError("Failed to construct a canonical molecule from the input.")

        self._mol = canonical_mol
        self.assign_atom_map(inplace=True, overwrite=overwrite_atom_map)
        self._formula = Formula.from_mol(self._mol)

    def __repr__(self) -> str:
        """
        Convert the compound to debug representation.

        Returns
        -------
        str
            Debug-style representation of the compound.
        """
        return f"Compound(smiles={self.smiles})"

    def __str__(self) -> str:
        """
        Convert the compound to SMILES.

        Returns
        -------
        str
            Canonical SMILES string without atom map numbers.
        """
        return self.smiles
    
    @classmethod
    def from_smiles(cls, smiles: str, overwrite_atom_map: bool = False) -> "Compound":
        """
        Create a compound from a SMILES string.

        Parameters
        ----------
        smiles : str
            Input SMILES string.
        overwrite_atom_map : bool, default=False
            Whether to overwrite existing atom map numbers.

        Returns
        -------
        Compound
            Compound created from the SMILES string.

        Raises
        ------
        ValueError
            If the SMILES string is invalid.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        return cls(mol, overwrite_atom_map=overwrite_atom_map)
    
    @property
    def mol(self) -> Chem.Mol:
        """
        RDKit molecule without atom map numbers.

        Returns
        -------
        Chem.Mol
            Copy of the internal molecule with atom map numbers cleared.
        """
        mol = Chem.Mol(self._mol)  # Create a copy to avoid modifying the original
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)  # Reset atom map numbers to 0
        return mol

    @property
    def mol_with_atom_map(self) -> Chem.Mol:
        """
        RDKit molecule with atom map numbers.

        Returns
        -------
        Chem.Mol
            Copy of the internal molecule including atom map numbers.
        """
        return Chem.Mol(self._mol)
    
    @property
    def smiles(self) -> str:
        """
        Canonical SMILES string without atom map numbers.

        Returns
        -------
        str
            Canonical SMILES string.
        """
        return Chem.MolToSmiles(self.mol, canonical=True)
    

    @property
    def mapped_smiles(self) -> str:
        """
        Canonical SMILES string with atom map numbers.

        Returns
        -------
        str
            Canonical SMILES string including atom map numbers.
        """
        return Chem.MolToSmiles(self._mol, canonical=True)
    
    @property
    def atom_index_mapped_smiles(self) -> str:
        """
        SMILES string with atom indices used as atom map numbers.

        Returns
        -------
        str
            SMILES string in which each atom map number equals its atom index.
        """
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise ValueError("Failed to reconstruct molecule from canonical SMILES.")

        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())

        return Chem.MolToSmiles(mol, canonical=True)
    
    @property
    def formula(self) -> Formula:
        """
        Molecular formula of the compound.

        Returns
        -------
        Formula
            Molecular formula.
        """
        return self._formula.copy()
    
    @property
    def atom_map_to_index(self) -> bidict[int, int]:
        """
        Mapping from atom map number to canonical atom index.

        Returns
        -------
        bidict[int, int]
            Bidirectional mapping from atom map number to atom index.

        Notes
        -----
        Atom map numbers are read from the internal molecule before
        canonicalization, then matched to the canonical atom order.
        """
        mol = Chem.Mol(self._mol)  # Create a copy to avoid modifying the original

        # Preserve original atom map numbers
        old_atom_map_num = {atom.GetIdx(): atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}

        # Remove atom map numbers for canonicalization
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

        # Generate canonical SMILES (this triggers atom ordering)
        Chem.MolToSmiles(mol, canonical=True)

        # Retrieve canonical atom order
        atom_order = [int(s) for s in mol.GetProp("_smilesAtomOutputOrder").replace('[', '').replace(']', '').split(",") if s !='']
        assert len(atom_order) == len(old_atom_map_num), "Atom count mismatch after canonicalization"

        # Build mapping: atom_map_num → new atom index
        atom_map_to_idx = {}
        for new_idx, old_idx in enumerate(atom_order):
            atom_map_num = old_atom_map_num.get(old_idx, None)
            if atom_map_num is not None:
                atom_map_to_idx[atom_map_num] = new_idx

        return bidict(atom_map_to_idx)
    
    @property
    def charge(self) -> int:
        """
        Formal charge of the compound.

        Returns
        -------
        int
            Sum of formal charges over all atoms.
        """
        return sum(atom.GetFormalCharge() for atom in self._mol.GetAtoms())
    
    @property
    def exact_mass(self) -> float:
        """
        Monoisotopic exact mass of the compound.

        Returns
        -------
        float
            Exact mass derived from the molecular formula.
        """
        return self.formula.exact_mass

    def assign_atom_map(
        self,
        inplace: bool = False,
        overwrite: bool = False,
        atom_map_dict: Optional[Dict[int, int]] = None,
    ) -> Optional["Compound"]:
        """
        Assign atom map numbers to atoms.

        Parameters
        ----------
        inplace : bool, default=False
            Whether to modify this compound in place.
        overwrite : bool, default=False
            Whether to overwrite existing atom map numbers.
        atom_map_dict : dict[int, int] or None, default=None
            Optional dictionary to store old-to-new atom map number mappings.

        Returns
        -------
        Compound or None
            Modified compound if ``inplace=False``, otherwise ``None``.

        Raises
        ------
        AssertionError
            If ``atom_map_dict`` is not empty when provided.
        """
        assert (atom_map_dict is None) or (isinstance(atom_map_dict, dict) and len(atom_map_dict) == 0), "atom_map_dict must be a dict or empty if provided"

        # Use a copy or the original mol
        mol = self._mol if inplace else Chem.Mol(self._mol)
        n_atoms = mol.GetNumAtoms()
        
        used_map_nums = set()
        if not overwrite:
            used_map_nums = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
            max_map_num = max(used_map_nums, default=0)
            diff = set(range(1, max_map_num)).difference(used_map_nums)
            next_map_nums = list(diff) + list(range(max_map_num + 1, max_map_num + n_atoms + 1))
        else:
            next_map_nums = list(range(1, n_atoms + 1))
        
        for atom in mol.GetAtoms():
            if overwrite or atom.GetAtomMapNum() == 0:
                old_map_num = atom.GetAtomMapNum()
                new_map_num = next_map_nums.pop(0)
                atom.SetAtomMapNum(new_map_num)
                if atom_map_dict is not None and old_map_num > 0:
                    atom_map_dict[old_map_num] = new_map_num

        if inplace:
            self._mol = mol
        else:
            return Compound(mol)

    def copy(self) -> "Compound":
        """
        Create a copy of the compound.

        Returns
        -------
        Compound
            Copied compound.
        """
        return Compound(self._mol)

    def get_atom_index_from_map(self, map_num: int) -> Optional[int]:
        """
        Get the atom index corresponding to an atom map number.

        Parameters
        ----------
        map_num : int
            Atom map number.

        Returns
        -------
        int or None
            Atom index if found, otherwise ``None``.
        """
        for atom in self._mol.GetAtoms():
            if atom.GetAtomMapNum() == map_num:
                return atom.GetIdx()
        return None
        


    