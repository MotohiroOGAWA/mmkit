from __future__ import annotations

import re
from collections import OrderedDict
from typing import Dict, Iterable, List, Union

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class Neutron:
    """Virtual representation of a neutron-like pseudo element."""

    symbol = "n"
    mass = 0.9987508525
    charge = 0
    atomic_num = 0


neutron = Neutron()


class Formula:
    """
    Chemical formula with element counts and net charge.

    Notes
    -----
    - Elements are stored in Hill order.
    - A virtual neutron symbol ("n") is supported.
    - Negative element counts are allowed for formula arithmetic.
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
        Create a Formula object.

        Parameters
        ----------
        elements:
            Mapping from element symbol to count.
        charge:
            Net charge of the formula.
        raw_formula:
            Original formula string, stored for reference.
        """
        self._charge: int = charge
        self._raw_formula: str = raw_formula
        self._elements: OrderedDict[str, int] = OrderedDict()
        self._reorder_elements(elements.copy())

    @classmethod
    def parse(cls, formula_str: str, store_raw: bool = False) -> "Formula":
        """
        Parse a formula string and return a Formula object.

        Parameters
        ----------
        formula_str:
            Formula string such as ``C6H12O6`` or ``C6H12O6+``.
        store_raw:
            Whether to store the original input in ``raw_formula``.

        Returns
        -------
        Formula
            Parsed Formula instance.
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
        Create a Formula from an RDKit Mol object.

        Implicit hydrogens are included by converting them to explicit
        hydrogens before counting atoms.
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
        """Return an empty formula."""
        return Formula(elements={}, charge=0)

    @staticmethod
    def neutron() -> Neutron:
        """Return the virtual neutron object."""
        return Neutron()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
    @property
    def elements(self) -> Dict[str, int]:
        """Return a copy of the element-count mapping."""
        return OrderedDict(self._elements)

    @property
    def charge(self) -> int:
        """Return the net charge."""
        return self._charge

    @property
    def raw_formula(self) -> str:
        """Return the original formula string."""
        return self._raw_formula

    @property
    def exact_mass(self) -> float:
        """
        Calculate the exact mass.

        Returns
        -------
        float
            Exact mass based on the most common isotope mass of each element.
        """
        periodic_table = Chem.GetPeriodicTable()
        mass = 0.0

        for elem, count in self._elements.items():
            if elem == neutron.symbol:
                mass += neutron.mass * count
            else:
                atomic_number = periodic_table.GetAtomicNumber(elem)
                mass += periodic_table.GetMostCommonIsotopeMass(atomic_number) * count

        return mass

    @property
    def is_nonnegative(self) -> bool:
        """
        Return whether all element counts are non-negative.
        """
        return all(count >= 0 for count in self._elements.values())

    @property
    def value(self) -> str:
        """Return the formula string including charge."""
        return self.to_string(no_charge=False)

    @property
    def plain_value(self) -> str:
        """Return the formula string without charge."""
        return self.to_string(no_charge=True)

    @property
    def plain(self) -> "Formula":
        """Return a copy without charge."""
        return Formula(self._elements, 0, self._raw_formula)

    @property
    def normalized(self) -> "Formula":
        """Return a copy with ``raw_formula`` cleared."""
        return Formula(self._elements.copy(), self._charge, "")

    @property
    def normalized_plain(self) -> "Formula":
        """Return a copy without charge and without ``raw_formula``."""
        return Formula(self._elements.copy(), 0, "")

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        return f"Formula({self})"

    def __str__(self) -> str:
        return self.to_string(no_charge=False)

    def to_string(self, no_charge: bool = False) -> str:
        """
        Convert the formula to a string.

        Parameters
        ----------
        no_charge:
            If True, omit the charge suffix.

        Returns
        -------
        str
            Formula string.
        """
        parts: List[str] = []
        neutron_part = ""

        for elem, count in self._elements.items():
            if elem == neutron.symbol:
                if count > 0:
                    neutron_part = f"+{elem}{count if count != 1 else ''}"
                elif count < 0:
                    neutron_part = f"-{elem}{-count if count != -1 else ''}"
                continue

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
        if not isinstance(other, Formula):
            return False
        return (
            self._charge == other._charge
            and self._elements == other._elements
        )

    def __hash__(self) -> int:
        return hash((tuple(self._elements.items()), self._charge, self._raw_formula))

    def eq(self, other: "Formula", ignore_raw: bool = True) -> bool:
        """
        Compare two Formula objects.

        Parameters
        ----------
        other:
            Formula to compare with.
        ignore_raw:
            If False, also compare ``raw_formula``.

        Returns
        -------
        bool
            True if equal.
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
        """Return the sum of two formulas."""
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
        """Return the difference of two formulas."""
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
        """Support reversed multiplication."""
        return self.__mul__(factor)

    # ------------------------------------------------------------------
    # Predicates / membership
    # ------------------------------------------------------------------
    def __contains__(self, item: Union["Formula", str]) -> bool:
        """
        Return whether this formula contains another formula.

        Only non-negative formulas are supported for containment checks.
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
        """Return a shallow copy of the Formula."""
        return Formula(self._elements, self._charge, self._raw_formula)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _reorder_elements(self, element_counts: Dict[str, int]) -> None:
        """
        Reorder elements according to the Hill system and remove zero counts.
        """
        ordered = self._reorder_element_keys(element_counts.keys())
        self._elements = OrderedDict(
            (k, element_counts[k]) for k in ordered if element_counts[k] != 0
        )

    @staticmethod
    def _split_formula_and_charge(formula: str) -> tuple[str, int]:
        """
        Split a formula string into body and charge.

        Examples
        --------
        ``C6H12O6+``  -> (``C6H12O6``, 1)
        ``Cl-``       -> (``Cl``, -1)
        ``Fe+2``      -> (``Fe``, 2)
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
        Parse element counts from a formula body without charge.
        """
        matches = re.findall(
            rf"([+-]?)([A-Z][a-z]?)(\d*)|([+-])({neutron.symbol})(\d*)",
            formula,
        )

        temp: Dict[str, int] = {}

        for sign, elem, count, n_sign, n_elem, n_count in matches:
            if n_elem == neutron.symbol:
                elem = neutron.symbol
                sign = n_sign
                count = n_count

            value = int(count) if count else 1
            if sign == "-":
                value = -value

            temp[elem] = temp.get(elem, 0) + value

        return {k: v for k, v in temp.items() if v != 0}

    @staticmethod
    def _reorder_element_keys(elements: Iterable[str]) -> tuple[str, ...]:
        """
        Reorder element symbols according to the Hill system.
        """
        mol = Chem.RWMol()
        exist_neutron = False

        for elem in elements:
            if elem == neutron.symbol:
                exist_neutron = True
                continue

            atom = Chem.Atom(elem)
            atom.SetNoImplicit(True)
            mol.AddAtom(atom)

        formula_str = rdMolDescriptors.CalcMolFormula(mol.GetMol())
        matches = re.findall(r"([A-Z][a-z]?)(\d*)", formula_str)

        ordered = tuple(match[0] for match in matches)
        if exist_neutron:
            ordered += (neutron.symbol,)

        return ordered