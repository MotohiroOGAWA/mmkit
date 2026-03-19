from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Literal, Tuple, Union, overload
import re

from rdkit import Chem  # noqa: F401

from .Formula import Formula
from ._parsing import charge_from_str


class Adduct:
    """
    Representation of an adduct ion.

    An ``Adduct`` describes how a neutral molecule is transformed into an
    observed ion by tracking:

    - ion type (e.g. ``"M"``, ``"F"``)
    - number of molecules
    - added formulas (adducts_in)
    - removed formulas (adducts_out)
    - net charge

    Examples
    --------
    ``[M+H]+``
        Protonated molecule.

    ``[2M+Na]+``
        Sodium-adducted dimer.

    ``[M-H]-``
        Deprotonated molecule.
    """

    def __init__(
        self,
        ion_type: Literal["M", "F"],
        n_molecules: int = 1,
        adducts_in: List[Formula] | None = None,
        adducts_out: List[Formula] | None = None,
        charge: int = 0,
    ):
        """
        Initialize an Adduct.

        Parameters
        ----------
        ion_type : {"M", "F"}
            Ion type identifier.
        n_molecules : int, default=1
            Number of neutral molecules.
        adducts_in : list[Formula] or None
            Formulas added during ionization.
        adducts_out : list[Formula] or None
            Formulas removed during ionization.
        charge : int, default=0
            Net charge of the ion.

        Notes
        -----
        Added and removed formulas are internally merged into a signed count
        representation.
        """
        assert ion_type in ["M", "F"], (
            f"ion_type must be one of 'M' or 'F', but got '{ion_type}'."
        )
        assert n_molecules >= 1, "n_molecules must be at least 1"

        adducts_in = [] if adducts_in is None else adducts_in
        adducts_out = [] if adducts_out is None else adducts_out

        self._ion_type = ion_type
        self._n_molecules = n_molecules

        formula_count: Dict[Formula, int] = defaultdict(int)

        for formula in adducts_in:
            formula_count[formula.plain] += 1

        for formula in adducts_out:
            formula_count[formula.plain] -= 1

        self._adduct_formulas: Dict[Formula, int] = {
            formula.copy(): count
            for formula, count in formula_count.items()
            if count != 0
        }
        self._charge = charge

    @classmethod
    def parse(cls, adduct_str: str) -> "Adduct":
        """
        Parse an adduct string.

        Parameters
        ----------
        adduct_str : str
            Adduct string such as ``"[M+H]+"`` or ``"[2M-H]-"``.

        Returns
        -------
        Adduct
            Parsed adduct object.

        Raises
        ------
        AssertionError
            If the string does not follow bracketed adduct notation.
        ValueError
            If the molecule identifier cannot be parsed.
        """
        assert adduct_str.startswith("[") and "]" in adduct_str, (
            f"Invalid adduct format: {adduct_str}"
        )

        main, charge_part = adduct_str[1:].split("]")
        charge = charge_from_str(charge_part.strip())

        n_match = re.match(r"(\d*)([A-Za-z]+)", main)
        if not n_match:
            raise ValueError(
                f"Invalid adduct format: missing molecule identifier in {adduct_str}"
            )

        n_molecules = int(n_match.group(1)) if n_match.group(1) else 1
        ion_type = n_match.group(2)
        remainder = main[n_match.end():]

        pattern = re.compile(r"([+-])(\d*)([A-Z][A-Za-z0-9]*)|([+-])(\d*)(i)")
        adducts_in: list[Formula] = []
        adducts_out: list[Formula] = []

        for sign, num, formula_str, i_sign, i_num, i_formula_str in pattern.findall(remainder):
            if i_formula_str == "i":
                formula_str = f"+{Formula.neutron().symbol}"
                sign = i_sign
                num = i_num

            count = int(num) if num else 1
            formulas = [Formula.parse(formula_str) for _ in range(count)]

            if sign == "+":
                adducts_in.extend(formulas)
            else:
                adducts_out.extend(formulas)

        adduct = cls(
            ion_type=ion_type,
            n_molecules=n_molecules,
            adducts_in=adducts_in,
            adducts_out=adducts_out,
            charge=charge,
        )
        return adduct
    
    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

    @property
    def ion_type(self) -> Literal["M", "F"]:
        """
        Ion type of the adduct.

        Returns
        -------
        Literal["M", "F"]
            Ion type of this adduct.
        """
        return self._ion_type

    @property
    def n_molecules(self) -> int:
        """
        Number of neutral molecules.

        Returns
        -------
        int
            Number of molecules represented by this adduct.
        """
        return self._n_molecules

    @property
    def charge(self) -> int:
        """
        Net charge of the adduct.

        Returns
        -------
        int
            Net charge of the adduct.
        """
        return self._charge

    @property
    def adduct_formulas(self) -> Dict[Formula, int]:
        """
        Signed counts of adduct formulas.

        Returns
        -------
        dict[Formula, int]
            Mapping from formula to signed count.
        """
        return dict(self._adduct_formulas)

    @property
    def formula_diff(self) -> Formula:
        """
        Net formula difference introduced by the adduct.

        Returns
        -------
        Formula
            Net added or removed composition.
        """
        total = Formula.empty()
        for formula, count in self._adduct_formulas.items():
            total = total + (formula * count)
        return total

    @property
    def mass_shift(self) -> float:
        """
        Exact mass shift introduced by the adduct.

        Returns
        -------
        float
            Total exact mass added or removed.
        """
        return sum(
            formula.exact_mass * count
            for formula, count in self._adduct_formulas.items()
        )

    @property
    def element_diff(self) -> dict[str, int]:
        """
        Total element count difference introduced by the adduct.

        Returns
        -------
        dict[str, int]
            Mapping from element symbol to signed count.
        """
        total: dict[str, int] = defaultdict(int)

        for formula, count in self._adduct_formulas.items():
            for elem, elem_count in formula.elements.items():
                total[elem] += elem_count * count

        return dict(total)

    # -------------------------------------------------------------------------
    # Basic methods
    # -------------------------------------------------------------------------

    def copy(self) -> "Adduct":
        """
        Create a copy of the adduct.

        Returns
        -------
        Adduct
            Copied adduct object.
        """
        adducts_in, adducts_out = self._split_adduct_formulas()

        copied = Adduct(
            ion_type=self._ion_type,
            n_molecules=self._n_molecules,
            adducts_in=[f.copy() for f in adducts_in],
            adducts_out=[f.copy() for f in adducts_out],
            charge=self._charge,
        )
        return copied

    def set_charge(self, charge: int) -> None:
        """
        Set the net charge of the adduct.

        Parameters
        ----------
        charge : int
            New charge value.
        """
        self._charge = charge

    @overload
    def get_formula_count(self, formula: Formula) -> int: ...

    @overload
    def get_formula_count(self, formula: str) -> int: ...

    def get_formula_count(self, formula: Union[Formula, str]) -> int:
        """
        Get the signed count of a specific formula.

        Parameters
        ----------
        formula : Formula or str
            Target formula.

        Returns
        -------
        int
            Positive if added, negative if removed, zero if absent.

        Raises
        ------
        TypeError
            If ``formula`` is neither ``Formula`` nor ``str``.
        """
        if isinstance(formula, str):
            formula = Formula.parse(formula)
        elif not isinstance(formula, Formula):
            raise TypeError(f"formula must be Formula or str, got {type(formula)}")

        return self._adduct_formulas.get(formula.plain, 0)

    # -------------------------------------------------------------------------
    # String representation
    # -------------------------------------------------------------------------

    def __str__(self) -> str:
        """
        Convert the adduct to canonical string form.

        Returns
        -------
        str
            Adduct notation such as ``[M+H]+`` or ``[M-H]-``.
        """
        parts: list[str] = []

        for formula, count in sorted(
            self._adduct_formulas.items(),
            key=lambda x: (-x[0].exact_mass, x[0].raw_formula),
        ):
            formula_str = formula.raw_formula.replace("+", "").replace("-", "")

            if count > 0:
                prefix = f"+{count}" if count > 1 else "+"
                parts.append(f"{prefix}{formula_str}")
            elif count < 0:
                prefix = f"-{abs(count)}" if count < -1 else "-"
                parts.append(f"{prefix}{formula_str}")

        body = "".join(parts)
        nM = f"{self._n_molecules if self._n_molecules > 1 else ''}{self._ion_type}"

        if self._charge > 0:
            charge = f"+{self._charge}" if self._charge > 1 else "+"
        elif self._charge < 0:
            charge = f"{self._charge}" if self._charge < -1 else "-"
        else:
            charge = ""

        return f"[{nM}{body}]{charge}"

    def __repr__(self) -> str:
        """
        Convert the adduct to debug representation.

        Returns
        -------
        str
            Debug-style representation of the adduct.
        """
        return f"Adduct({self})"

    def __eq__(self, other: object) -> bool:
        """
        Compare two adducts for equality.

        Parameters
        ----------
        other : object
            Object to compare.

        Returns
        -------
        bool
            True if the canonical string representations are equal.
        """
        return str(self) == str(other)

    def __hash__(self) -> int:
        """
        Compute the hash value of the adduct.

        Returns
        -------
        int
            Hash based on the canonical string representation.
        """
        return hash(str(self))

    # -------------------------------------------------------------------------
    # Adduct combination
    # -------------------------------------------------------------------------

    def add(
        self,
        other: "Adduct",
        prefer_ion_type: bool = False,
        prefer_n_molecules: bool = False,
        prefer_charge: bool = False,
    ) -> "Adduct":
        """
        Combine this adduct with another adduct.

        Parameters
        ----------
        other : Adduct
            Adduct to combine with.
        prefer_ion_type : bool, default=False
            Whether to keep ``self.ion_type`` when ion types differ.
        prefer_n_molecules : bool, default=False
            Whether to keep ``self.n_molecules`` when molecule counts differ.
        prefer_charge : bool, default=False
            Whether to keep ``self.charge`` when charges differ.

        Returns
        -------
        Adduct
            Combined adduct.

        Raises
        ------
        AssertionError
            If ``other`` is not an ``Adduct``.
        ValueError
            If ion types differ and ``prefer_ion_type`` is False.
        """
        assert isinstance(other, Adduct), (
            f"Can only add Adduct to Adduct, but got {type(other)}"
        )

        if not prefer_ion_type and self._ion_type != other._ion_type:
            raise ValueError(
                f"Cannot add Adducts with different ion_type: "
                f"{self._ion_type} + {other._ion_type}"
            )

        if not prefer_n_molecules:
            assert self._n_molecules == other._n_molecules, (
                f"Cannot add Adducts with different n_molecules: "
                f"{self._n_molecules} + {other._n_molecules}"
            )

        if not prefer_charge:
            assert self._charge == other._charge, (
                f"Cannot add Adducts with different charge: "
                f"{self._charge} + {other._charge}"
            )

        merged: dict[Formula, int] = defaultdict(int)

        for formula, count in self._adduct_formulas.items():
            merged[formula] += count
        for formula, count in other._adduct_formulas.items():
            merged[formula] += count

        adducts_in: list[Formula] = []
        adducts_out: list[Formula] = []

        for formula, count in merged.items():
            if count > 0:
                adducts_in.extend([formula.copy()] * count)
            elif count < 0:
                adducts_out.extend([formula.copy()] * (-count))

        out = Adduct(
            ion_type=self._ion_type,
            n_molecules=self._n_molecules,
            adducts_in=adducts_in,
            adducts_out=adducts_out,
            charge=self._charge,
        )

        return out

    def add_prefer_self(self, other: "Adduct") -> "Adduct":
        """
        Combine two adducts while keeping this adduct's metadata.

        Parameters
        ----------
        other : Adduct
            Adduct to combine with.

        Returns
        -------
        Adduct
            Combined adduct.
        """
        return self.add(
            other,
            prefer_ion_type=True,
            prefer_n_molecules=True,
            prefer_charge=True,
        )

    # -------------------------------------------------------------------------
    # Formula and mass calculations
    # -------------------------------------------------------------------------

    def apply_to_formula(self, neutral_formula: Formula) -> Formula:
        """
        Apply the adduct to a neutral formula.

        Parameters
        ----------
        neutral_formula : Formula
            Neutral molecular formula.

        Returns
        -------
        Formula
            Resulting ion formula.
        """
        total_formula = neutral_formula * self._n_molecules + self.formula_diff
        total_formula._charge = self.charge
        return total_formula

    def apply_to_mass(self, neutral_mass: float) -> float:
        """
        Apply the adduct to a neutral mass.

        Parameters
        ----------
        neutral_mass : float
            Neutral monoisotopic mass.

        Returns
        -------
        float
            Ion mass before charge division.
        """
        return neutral_mass * self._n_molecules + self.mass_shift

    def apply_to_mz(self, neutral_mass: float) -> float:
        """
        Apply the adduct to a neutral mass and compute m/z.

        Parameters
        ----------
        neutral_mass : float
            Neutral monoisotopic mass.

        Returns
        -------
        float
            Observed m/z value.

        Raises
        ------
        ValueError
            If the adduct charge is zero.
        """
        if self.charge == 0:
            raise ValueError(f"Cannot calculate m/z for uncharged adduct: {self}")

        return self.apply_to_mass(neutral_mass) / abs(self.charge)

    # -------------------------------------------------------------------------
    # Splitting helpers
    # -------------------------------------------------------------------------

    def _split_adduct_formulas(self) -> Tuple[List[Formula], List[Formula]]:
        """
        Split signed formula counts into input and output lists.

        Returns
        -------
        Tuple[List[Formula], List[Formula]]
            Tuple of ``(adducts_in, adducts_out)``.
        """
        adducts_in: List[Formula] = []
        adducts_out: List[Formula] = []

        for formula, count in self._adduct_formulas.items():
            if count > 0:
                adducts_in.extend([formula.copy()] * count)
            elif count < 0:
                adducts_out.extend([formula.copy()] * (-count))

        return adducts_in, adducts_out

    def split(self, split_each: bool = False) -> Tuple["Adduct", ...]:
        """
        Split the adduct into independent components.

        Parameters
        ----------
        split_each : bool, default=False
            Whether to split each formula into its own adduct.

        Returns
        -------
        Tuple[Adduct, ...]
            Split adduct objects.
        """
        result: List[Adduct] = []
        adducts_in, adducts_out = self._split_adduct_formulas()

        if adducts_in:
            if split_each:
                for formula in adducts_in:
                    result.append(
                        Adduct(
                            ion_type=self._ion_type,
                            n_molecules=self._n_molecules,
                            adducts_in=[formula.copy()],
                            adducts_out=[],
                            charge=0,
                        )
                    )
            else:
                result.append(
                    Adduct(
                        ion_type=self._ion_type,
                        n_molecules=self._n_molecules,
                        adducts_in=adducts_in,
                        adducts_out=[],
                        charge=0,
                    )
                )

        if adducts_out:
            if split_each:
                for formula in adducts_out:
                    result.append(
                        Adduct(
                            ion_type=self._ion_type,
                            n_molecules=self._n_molecules,
                            adducts_in=[],
                            adducts_out=[formula.copy()],
                            charge=0,
                        )
                    )
            else:
                result.append(
                    Adduct(
                        ion_type=self._ion_type,
                        n_molecules=self._n_molecules,
                        adducts_in=[],
                        adducts_out=adducts_out,
                        charge=0,
                    )
                )

        return tuple(result)

    @staticmethod
    def decompose_adduct(
        adduct: "Adduct",
        reference_adducts: Tuple["Adduct", ...],
    ) -> Tuple[Dict["Adduct", int], "Adduct"]:
        """
        Decompose an adduct into known adduct types and a neutral component.

        Parameters
        ----------
        adduct : Adduct
            Target adduct to decompose.
        reference_adducts : tuple[Adduct, ...]
            Supported adduct types used for matching.

        Returns
        -------
        tuple[dict[Adduct, int], Adduct]
            - Mapping from supported adduct type to count.
            - Remaining neutral component (charge = 0).

        Raises
        ------
        ValueError
            If an unsupported adduct formula is encountered.
        NotImplementedError
            If negative adduct decomposition is requested.
        """
        from collections import defaultdict

        # Frequently used neutral formulas
        h2o = Formula.parse("H2O")
        neutron = Formula.parse("+n")

        # Result: supported adduct composition
        adduct_composition: Dict[Adduct, int] = defaultdict(int)

        # Neutral components (added / removed)
        neutral_in: List[Formula] = []
        neutral_out: List[Formula] = []

        # ------------------------------------------------------------------
        # Positive mode
        # ------------------------------------------------------------------
        if adduct.charge > 0:
            # Build lookup: formula string → supported adduct
            supported_lookup = {
                str(ref.formula_diff): ref
                for ref in reference_adducts
                if ref.charge > 0
            }

            # Iterate over all sub-formulas in the adduct
            for formula, count in adduct._adduct_formulas.items():
                if count == 0:
                    continue

                formula_str = str(formula)

                # ----------------------------------------------------------
                # Case 1: supported adduct type
                # ----------------------------------------------------------
                if formula_str in supported_lookup:
                    adduct_type = supported_lookup[formula_str]
                    adduct_composition[adduct_type] += count
                    continue

                # ----------------------------------------------------------
                # Case 2: neutral or special species
                # ----------------------------------------------------------
                if count > 0:
                    if formula in (h2o, neutron):
                        neutral_in.extend([formula.copy()] * count)
                    else:
                        raise ValueError(
                            f"Unsupported positive adduct formula: {formula} in {adduct}"
                        )
                else:
                    neutral_out.extend([formula.copy()] * (-count))

        # ------------------------------------------------------------------
        # Negative mode (not implemented)
        # ------------------------------------------------------------------
        elif adduct.charge < 0:
            raise NotImplementedError("Negative adduct decomposition is not implemented.")

        # ------------------------------------------------------------------
        # Neutral (invalid)
        # ------------------------------------------------------------------
        else:
            raise ValueError(f"Adduct must have non-zero charge: {adduct}")

        # ------------------------------------------------------------------
        # Build neutral component (charge = 0)
        # ------------------------------------------------------------------
        neutral_component = Adduct(
            ion_type=adduct._ion_type,
            n_molecules=adduct._n_molecules,
            adducts_in=neutral_in,
            adducts_out=neutral_out,
            charge=0,
        )

        # Sort output for deterministic behavior
        adduct_composition = dict(
            sorted(adduct_composition.items(), key=lambda x: str(x[0]))
        )

        return adduct_composition, neutral_component