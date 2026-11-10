using System;
using Mmkit;

using var compound = Compound.from_smiles("CCO");
using var formula = compound.formula();
using var adduct = Adduct.parse("[M+H]+");
if (formula.to_string() != "C2H6O") throw new Exception("formula");
if (Math.Abs(adduct.apply_to_mz(compound.exact_mass()) - 47.049689812) > 1e-6)
    throw new Exception("m/z");
bool rejected = false;
try { using var invalid = Adduct.parse("bad"); }
catch (ArgumentException) { rejected = true; }
if (!rejected) throw new Exception("invalid input accepted");
Console.WriteLine("C# P/Invoke OK: " + formula.to_string());
