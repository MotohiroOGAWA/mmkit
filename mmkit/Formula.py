from __future__ import annotations

import re
from collections import OrderedDict
from typing import Dict, Iterable, List, Union

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class Formula:
    """
    Represent a chemical formula with element counts and net charge.

    A ``Formula`` stores an elemental composition together with a formal charge.
    It supports parsing from strings, arithmetic operations between formulas,
    exact-mass calculation, and conversion back to canonical string form.

    Notes
    -----
    - Elements are stored in Hill order.
    - Negative element counts are allowed for formula arithmetic.
    - ``raw_formula`` can preserve the original input string when needed.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    def __init__(
        self,
        elements: Dict[str, int],
        charge: int = 0,
        raw_formula: str = "",
    ) -> None:
        """
        Initialize a Formula.

        Parameters
        ----------
        elements
            Mapping from element symbol to count.
        charge
            Net charge of the formula.
        raw_formula
            Original formula string stored for reference.

        Notes
        -----
        Element keys are reordered into Hill order and zero-count elements
        are removed during initialization.
        """
        self._charge: int = charge
        self._raw_formula: str = raw_formula
        self._elements: OrderedDict[str, int] = OrderedDict()
        self._reorder_elements(elements.copy())

    @classmethod
    def parse(cls, formula_str: str, store_raw: bool = False) -> "Formula":
        """
        Parse a formula string into a Formula object.

        Parameters
        ----------
        formula_str
            Formula string such as ``C6H12O6`` or ``C6H12O6+``.
        store_raw
            Whether to store the original input string in ``raw_formula``.

        Returns
        -------
        Formula
            Parsed Formula instance.

        Examples
        --------
        >>> Formula.parse("HOH").value
        'H2O'
        >>> Formula.parse("C6H12O6+").charge
        1
        """
        formula_body, charge = cls._split_formula_and_charge(formula_str)
        elements = cls._parse_element_counts(formula_body)

        return cls(
            elements=elements,
            charge=charge,
            raw_formula=formula_str if store_raw else "",
        )

    @classmethod
    def from_mol(cls, mol: Chem.Mol) -> "Formula":
        """
        Create a Formula from an RDKit molecule.

        Parameters
        ----------
        mol
            RDKit molecule object.

        Returns
        -------
        Formula
            Formula derived from the molecule.

        Notes
        -----
        Implicit hydrogens are converted to explicit hydrogens before counting
        atoms, so hydrogen counts are included in the result.
        """
        mol_with_h = Chem.AddHs(mol)

        elements: Dict[str, int] = {}
        total_charge = 0

        for atom in mol_with_h.GetAtoms():
            symbol = atom.GetSymbol()
            elements[symbol] = elements.get(symbol, 0) + 1
            total_charge += atom.GetFormalCharge()

        return cls(elements=elements, charge=total_charge)

    @staticmethod
    def empty() -> "Formula":
        """
        Return an empty formula.

        Returns
        -------
        Formula
            Formula with no elements and zero charge.
        """
        return Formula(elements={}, charge=0)

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def elements(self) -> Dict[str, int]:
        """
        Element-count mapping in Hill order.

        Returns
        -------
        OrderedDict[str, int]
            Copy of the internal mapping from element symbol to count.
        """
        return OrderedDict(self._elements)

    @property
    def charge(self) -> int:
        """
        Net charge of the formula.

        Returns
        -------
        int
            Formal charge of the formula.
        """
        return self._charge

    @property
    def raw_formula(self) -> str:
        """
        Original formula string, if available.

        Returns
        -------
        str
            Stored original formula string if available. If no raw string was
            stored, the canonical formula string is returned instead.
        """
        if self._raw_formula == "":
            return self.to_string(no_charge=False)
        return self._raw_formula

    @property
    def exact_mass(self) -> float:
        """
        Monoisotopic exact mass of the formula.

        Returns
        -------
        float
            Exact mass calculated from the most common isotope mass of each
            element.
        """
        periodic_table = Chem.GetPeriodicTable()
        mass = 0.0

        for elem, count in self._elements.items():
            atomic_number = periodic_table.GetAtomicNumber(elem)
            mass += periodic_table.GetMostCommonIsotopeMass(atomic_number) * count

        return mass

    @property
    def is_nonnegative(self) -> bool:
        """
        Whether all element counts are non-negative.

        Returns
        -------
        bool
            True if every element count is greater than or equal to zero.
        """
        return all(count >= 0 for count in self._elements.values())

    @property
    def value(self) -> str:
        """
        Canonical formula string including charge.

        Returns
        -------
        str
            Formula string with charge suffix.
        """
        return self.to_string(no_charge=False)

    @property
    def plain_value(self) -> str:
        """
        Canonical formula string without charge.

        Returns
        -------
        str
            Formula string without charge suffix.
        """
        return self.to_string(no_charge=True)

    @property
    def plain(self) -> "Formula":
        """
        Copy of the formula with zero charge.

        Returns
        -------
        Formula
            Copy of the formula with zero charge.
        """
        return Formula(self._elements, 0, self._raw_formula)

    @property
    def normalized(self) -> "Formula":
        """
        Copy with ``raw_formula`` cleared.

        Returns
        -------
        Formula
            Copy of the formula with the same composition and charge, but with
            an empty ``raw_formula``.
        """
        return Formula(self._elements.copy(), self._charge, "")

    @property
    def normalized_plain(self) -> "Formula":
        """
        Copy with zero charge and empty ``raw_formula``.

        Returns
        -------
        Formula
            Copy of the formula with zero charge and an empty ``raw_formula``.
        """
        return Formula(self._elements.copy(), 0, "")

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Return the debug representation.

        Returns
        -------
        str
            Debug-style representation of the formula.
        """
        return f"Formula({self})"

    def __str__(self) -> str:
        """
        Return the canonical formula string.

        Returns
        -------
        str
            Formula string including charge.
        """
        return self.to_string(no_charge=False)

    def to_string(self, no_charge: bool = False) -> str:
        """
        Convert the formula to a string.

        Parameters
        ----------
        no_charge
            If True, omit the charge suffix.

        Returns
        -------
        str
            Canonical formula string.

        Notes
        -----
        Elements are emitted in Hill order. Negative element counts are written
        with a leading ``-``.
        """
        parts: List[str] = []
        neutron_part = ""

        for elem, count in self._elements.items():
            if count > 0:
                parts.append(f"{elem}{count if count != 1 else ''}")
            elif count < 0:
                parts.append(f"-{elem}{-count if count != -1 else ''}")

        formula = "".join(parts) + neutron_part

        if not no_charge:
            if self._charge > 0:
                formula += "+" if self._charge == 1 else f"+{self._charge}"
            elif self._charge < 0:
                formula += "-" if self._charge == -1 else f"-{-self._charge}"

        return formula

    # ------------------------------------------------------------------
    # Comparison / hashing
    # ------------------------------------------------------------------
    def __eq__(self, other: object) -> bool:
        """
        Return whether two Formula objects are equal.

        Parameters
        ----------
        other
            Object to compare against.

        Returns
        -------
        bool
            True if both charge and element counts are equal.
        """
        if not isinstance(other, Formula):
            return False
        return (
            self._charge == other._charge
            and self._elements == other._elements
        )

    def __hash__(self) -> int:
        """
        Return the hash value.

        Returns
        -------
        int
            Hash based on element counts, charge, and raw formula string.
        """
        return hash((tuple(self._elements.items()), self._charge, self._raw_formula))

    def eq(self, other: "Formula", ignore_raw: bool = True) -> bool:
        """
        Compare two Formula objects with optional raw-string comparison.

        Parameters
        ----------
        other
            Formula to compare with.
        ignore_raw
            If False, also compare ``raw_formula``.

        Returns
        -------
        bool
            True if the two Formula objects are considered equal.
        """
        if not isinstance(other, Formula):
            return False

        base_equal = self == other
        if ignore_raw:
            return base_equal
        return base_equal and (self._raw_formula == other._raw_formula)

    # ------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------
    def __add__(self, other: "Formula") -> "Formula":
        """
        Return the sum of two formulas.

        Parameters
        ----------
        other
            Formula to add.

        Returns
        -------
        Formula
            Sum of the two formulas.
        """
        if not isinstance(other, Formula):
            return NotImplemented

        combined = dict(self._elements)
        for elem, count in other._elements.items():
            combined[elem] = combined.get(elem, 0) + count

        return Formula(
            elements=combined,
            charge=self._charge + other._charge,
        )

    def __sub__(self, other: "Formula") -> "Formula":
        """
        Return the difference of two formulas.

        Parameters
        ----------
        other
            Formula to subtract.

        Returns
        -------
        Formula
            Difference of the two formulas.
        """
        if not isinstance(other, Formula):
            return NotImplemented

        combined = dict(self._elements)
        for elem, count in other._elements.items():
            combined[elem] = combined.get(elem, 0) - count

        return Formula(
            elements=combined,
            charge=self._charge - other._charge,
        )

    def __mul__(self, factor: int) -> "Formula":
        """
        Multiply the formula by an integer.

        Parameters
        ----------
        factor
            Integer scaling factor.

        Returns
        -------
        Formula
            Scaled formula.

        Raises
        ------
        TypeError
            If ``factor`` is not an integer.

        Examples
        --------
        ``H2O * 2 -> H4O2``
        """
        if not isinstance(factor, int):
            raise TypeError(f"Formula can only be multiplied by int, not {type(factor)}")

        new_elements = {elem: count * factor for elem, count in self._elements.items()}
        new_charge = self._charge * factor
        return Formula(new_elements, new_charge, self._raw_formula)

    def __rmul__(self, factor: int) -> "Formula":
        """
        Support reversed multiplication.

        Parameters
        ----------
        factor
            Integer scaling factor.

        Returns
        -------
        Formula
            Scaled formula.
        """
        return self.__mul__(factor)

    # ------------------------------------------------------------------
    # Predicates / membership
    # ------------------------------------------------------------------
    def __contains__(self, item: Union["Formula", str]) -> bool:
        """
        Return whether this formula contains another formula.

        Parameters
        ----------
        item
            Formula object or formula string to test for containment.

        Returns
        -------
        bool
            True if this formula contains at least the counts required by
            ``item``.

        Raises
        ------
        ValueError
            If this formula contains negative element counts.
        TypeError
            If ``item`` is neither ``Formula`` nor ``str``.

        Notes
        -----
        Containment checks are only supported for non-negative formulas.
        """
        if not self.is_nonnegative:
            raise ValueError("Containment check is only supported for non-negative formulas.")

        if isinstance(item, str):
            item = Formula.parse(item)
        elif not isinstance(item, Formula):
            raise TypeError(f"Containment check only supports Formula or str, not {type(item)}")

        return all(self._elements.get(elem, 0) >= count for elem, count in item._elements.items())

    # ------------------------------------------------------------------
    # Copy utilities
    # ------------------------------------------------------------------
    def copy(self) -> "Formula":
        """
        Create a copy of the formula.

        Returns
        -------
        Formula
            Shallow copy of the Formula object.
        """
        return Formula(self._elements, self._charge, self._raw_formula)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _reorder_elements(self, element_counts: Dict[str, int]) -> None:
        """
        Reorder elements according to the Hill system.

        Parameters
        ----------
        element_counts
            Mapping from element symbol to count.

        Notes
        -----
        Zero-count elements are removed.
        """
        ordered = self._reorder_element_keys(element_counts.keys())
        self._elements = OrderedDict(
            (k, element_counts[k]) for k in ordered if element_counts[k] != 0
        )

    @staticmethod
    def _split_formula_and_charge(formula: str) -> tuple[str, int]:
        """
        Split a formula string into formula body and charge.

        Parameters
        ----------
        formula
            Formula string including an optional charge suffix.

        Returns
        -------
        tuple[str, int]
            Tuple of ``(formula_body, charge)``.

        Examples
        --------
        ``C6H12O6+`` -> (``C6H12O6``, 1)
        ``Cl-`` -> (``Cl``, -1)
        ``Fe+2`` -> (``Fe``, 2)
        """
        charge = 0
        charge_match = re.search(r"([+-]+|[+-]\d+)$", formula)

        if charge_match:
            charge_str = charge_match.group(1)
            formula = formula[: -len(charge_str)]

            charge = int(charge_str[1:]) if charge_str[1:] else 1
            if charge_str[0] == "-":
                charge *= -1

        return formula, charge

    @staticmethod
    def _parse_element_counts(formula: str) -> Dict[str, int]:
        """
        Parse element counts from a neutral formula body.

        Parameters
        ----------
        formula
            Formula string without a charge suffix.

        Returns
        -------
        dict[str, int]
            Mapping from element symbol to count.

        Notes
        -----
        Negative terms such as ``-H2`` are supported.
        """
        matches = re.findall(
            rf"([+-]?)([A-Z][a-z]?)(\d*)",
            formula,
        )

        temp: Dict[str, int] = {}

        for sign, elem, count in matches:
            value = int(count) if count else 1
            if sign == "-":
                value = -value

            temp[elem] = temp.get(elem, 0) + value

        return {k: v for k, v in temp.items() if v != 0}

    @staticmethod
    def _reorder_element_keys(elements: Iterable[str]) -> tuple[str, ...]:
        """
        Reorder element symbols according to the Hill system.

        Parameters
        ----------
        elements
            Iterable of element symbols.

        Returns
        -------
        tuple[str, ...]
            Element symbols in Hill order.
        """
        mol = Chem.RWMol()

        for elem in elements:

            atom = Chem.Atom(elem)
            atom.SetNoImplicit(True)
            mol.AddAtom(atom)

        formula_str = rdMolDescriptors.CalcMolFormula(mol.GetMol())
        matches = re.findall(r"([A-Z][a-z]?)(\d*)", formula_str)

        ordered = tuple(match[0] for match in matches)

        return ordered