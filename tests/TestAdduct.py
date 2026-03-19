import unittest

from mmkit.Formula import Formula
from mmkit.Adduct import Adduct


class TestAdduct(unittest.TestCase):
    """Unit tests for Adduct."""

    def test_init_basic(self):
        adduct = Adduct(
            ion_type="M",
            n_molecules=1,
            adducts_in=[Formula.parse("H+")],
            adducts_out=[],
            charge=1,
        )
        self.assertEqual(adduct.ion_type, "M")
        self.assertEqual(adduct.n_molecules, 1)
        self.assertEqual(adduct.charge, 1)

    def test_str_positive(self):
        adduct = Adduct(
            ion_type="M",
            adducts_in=[Formula.parse("H+")],
            charge=1,
        )
        self.assertEqual(str(adduct), "[M+H]+")

    def test_str_negative(self):
        adduct = Adduct(
            ion_type="M",
            adducts_out=[Formula.parse("H+")],
            charge=-1,
        )
        self.assertEqual(str(adduct), "[M-H]-")

    def test_str_multiple_molecules(self):
        adduct = Adduct(
            ion_type="M",
            n_molecules=2,
            adducts_in=[Formula.parse("Na+")],
            charge=1,
        )
        self.assertEqual(str(adduct), "[2M+Na]+")

    def test_repr(self):
        adduct = Adduct(
            ion_type="M",
            adducts_in=[Formula.parse("H+")],
            charge=1,
        )
        self.assertEqual(repr(adduct), "Adduct([M+H]+)")

    def test_eq_and_hash(self):
        adduct1 = Adduct.parse("[M+H]+")
        adduct2 = Adduct(
            ion_type="M",
            adducts_in=[Formula.parse("H+")],
            charge=1,
        )
        self.assertEqual(adduct1, adduct2)
        self.assertEqual(hash(adduct1), hash(adduct2))

    def test_adduct_formulas_returns_copy(self):
        adduct = Adduct.parse("[M+H]+")
        d = adduct.adduct_formulas
        d[Formula.parse("Na+")] = 100
        self.assertEqual(adduct.get_formula_count("Na+"), 0)

    def test_formula_diff(self):
        adduct = Adduct.parse("[M+H]+")
        self.assertEqual(adduct.formula_diff, Formula.parse("H"))

    def test_mass_shift(self):
        adduct = Adduct.parse("[M+H]+")
        self.assertAlmostEqual(
            adduct.mass_shift,
            Formula.parse("H+").exact_mass,
            places=6,
        )

    def test_element_diff(self):
        adduct = Adduct.parse("[M+H]+")
        self.assertEqual(adduct.element_diff["H"], 1)

        adduct2 = Adduct.parse("[M-H]-")
        self.assertEqual(adduct2.element_diff["H"], -1)

    def test_get_formula_count_with_formula(self):
        adduct = Adduct.parse("[M+H]+")
        self.assertEqual(adduct.get_formula_count(Formula.parse("H+")), 1)
        self.assertEqual(adduct.get_formula_count(Formula.parse("Na+")), 0)

    def test_get_formula_count_with_string(self):
        adduct = Adduct.parse("[M+H-H2O]+")
        self.assertEqual(adduct.get_formula_count("H+"), 1)
        self.assertEqual(adduct.get_formula_count("H2O"), -1)
        self.assertEqual(adduct.get_formula_count("Na+"), 0)

    def test_copy(self):
        adduct1 = Adduct.parse("[2M+Na]+")
        adduct2 = adduct1.copy()

        self.assertEqual(adduct1, adduct2)
        self.assertIsNot(adduct1, adduct2)
        self.assertEqual(adduct1.charge, adduct2.charge)

    def test_parse_positive(self):
        adduct = Adduct.parse("[M+H]+")
        self.assertEqual(adduct.ion_type, "M")
        self.assertEqual(adduct.n_molecules, 1)
        self.assertEqual(adduct.charge, 1)
        self.assertEqual(adduct.get_formula_count("H"), 1)

    def test_parse_negative(self):
        adduct = Adduct.parse("[M-H]-")
        self.assertEqual(adduct.ion_type, "M")
        self.assertEqual(adduct.n_molecules, 1)
        self.assertEqual(adduct.charge, -1)
        self.assertEqual(adduct.get_formula_count("H+"), -1)

    def test_parse_multiple_molecules(self):
        adduct = Adduct.parse("[2M+Na]+")
        self.assertEqual(adduct.ion_type, "M")
        self.assertEqual(adduct.n_molecules, 2)
        self.assertEqual(adduct.charge, 1)
        self.assertEqual(adduct.get_formula_count("Na+"), 1)

    def test_parse_complex(self):
        adduct = Adduct.parse("[M+HCOOH-H]-")
        self.assertEqual(adduct.charge, -1)
        self.assertEqual(adduct.get_formula_count("HCOOH"), 1)
        self.assertEqual(adduct.get_formula_count("H+"), -1)

    def test_parse_invalid_format(self):
        invalid_cases = [
            "M+H+",
            "",
            "abc",
            "[+H]",
        ]
        for s in invalid_cases:
            with self.subTest(adduct_str=s):
                with self.assertRaises((AssertionError, ValueError)):
                    Adduct.parse(s)

    def test_add_same_adduct_type(self):
        adduct1 = Adduct.parse("[M+H]+")
        adduct2 = Adduct.parse("[M+Na]+")
        adduct3 = adduct1.add(adduct2, prefer_charge=True)

        self.assertEqual(adduct3.ion_type, "M")
        self.assertEqual(adduct3.n_molecules, 1)
        self.assertEqual(adduct3.charge, 1)
        self.assertEqual(adduct3.get_formula_count("H"), 1)
        self.assertEqual(adduct3.get_formula_count("Na"), 1)

    def test_add_different_ion_type_raises(self):
        adduct1 = Adduct(ion_type="M")
        adduct2 = Adduct(ion_type="F")

        with self.assertRaises(ValueError):
            adduct1.add(adduct2)

    def test_add_different_n_molecules_raises(self):
        adduct1 = Adduct.parse("[M+H]+")
        adduct2 = Adduct.parse("[2M+H]+")

        with self.assertRaises(AssertionError):
            adduct1.add(adduct2)

    def test_add_different_charge_raises(self):
        adduct1 = Adduct.parse("[M+H]+")
        adduct2 = Adduct.parse("[M-H]-")

        with self.assertRaises(AssertionError):
            adduct1.add(adduct2)

    def test_add_prefer_self(self):
        adduct1 = Adduct.parse("[M+H]+")
        adduct2 = Adduct.parse("[2M+Na]+")

        adduct3 = adduct1.add_prefer_self(adduct2)
        self.assertEqual(adduct3.ion_type, "M")
        self.assertEqual(adduct3.n_molecules, 1)
        self.assertEqual(adduct3.charge, 1)
        self.assertEqual(adduct3.get_formula_count("H"), 1)
        self.assertEqual(adduct3.get_formula_count("Na"), 1)

    def test_apply_to_formula(self):
        adduct = Adduct.parse("[M+H]+")
        neutral = Formula.parse("C6H12O6")
        ion = adduct.apply_to_formula(neutral)

        expected = Formula.parse("C6H13O6+")
        self.assertEqual(ion, expected)
        self.assertEqual(ion.charge, 1)

    def test_apply_to_mass(self):
        adduct = Adduct.parse("[M+H]+")
        neutral_mass = 100.0
        self.assertAlmostEqual(
            adduct.apply_to_mass(neutral_mass),
            100.0 + Formula.parse("H+").exact_mass,
            places=6,
        )

    def test_apply_to_mz(self):
        adduct = Adduct.parse("[M+H]+")
        neutral_mass = 100.0
        expected = (100.0 + Formula.parse("H+").exact_mass) / 1
        self.assertAlmostEqual(adduct.apply_to_mz(neutral_mass), expected, places=6)

    def test_apply_to_mz_raises_for_zero_charge(self):
        adduct = Adduct(
            ion_type="M",
            adducts_in=[Formula.parse("H+")],
        )
        adduct.set_charge(0)

        with self.assertRaises(ValueError):
            adduct.apply_to_mz(100.0)

    def test_split_grouped(self):
        adduct = Adduct.parse("[M+H-H2O]+")
        parts = adduct.split(split_each=False)

        self.assertEqual(len(parts), 2)
        self.assertEqual(str(parts[0]), "[M+H]")
        self.assertEqual(str(parts[1]), "[M-H2O]")

    def test_split_each(self):
        adduct = Adduct.parse("[M+H+Na-H2O]+")
        parts = adduct.split(split_each=True)

        self.assertEqual(len(parts), 3)
        self.assertIn("[M+H]", [str(x) for x in parts])
        self.assertIn("[M+Na]", [str(x) for x in parts])
        self.assertIn("[M-H2O]", [str(x) for x in parts])


if __name__ == "__main__":
    unittest.main()