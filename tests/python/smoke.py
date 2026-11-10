from mmkit import Formula, Compound, Adduct

compound = Compound.from_smiles("CCO")
assert compound.formula == Formula.parse("C2H6O")
adduct = Adduct.parse("[M+H]+")
assert adduct.apply_to_formula(compound.formula).value == "C2H7O+"
assert abs(adduct.apply_to_mz(compound.exact_mass) - 47.049689812) < 1e-6
assert Formula.from_mol(compound.mol) == compound.formula
try:
    Adduct.parse("bad")
except ValueError:
    pass
else:
    raise AssertionError("Invalid input accepted")
print("Python bindings OK")
