#include "mmkit/Compound.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
template<class Action> void rejects(Action action) {
    try { action(); } catch (const std::exception&) { return; }
    throw std::runtime_error("Expected exception");
}
int main() {
    using mmkit::Compound;
    try {
        const auto ethanol = Compound::from_smiles("OCC");
        require(ethanol.smiles() == "CCO", "Canonical SMILES");
        require(ethanol.formula().to_string() == "C2H6O", "Formula");
        require(std::abs(ethanol.exact_mass() - 46.041865) < 0.000001, "Mass");
        require(ethanol.charge() == 0, "Neutral charge");
        require(ethanol.equals(Compound("CCO")), "Equivalent structures");
        require(!ethanol.equals(Compound("COC")), "Different structures");
        require(Compound("[NH4+]").charge() == 1, "Positive charge");
        require(Compound("[Cl-]").formula().to_string() == "Cl-", "Negative charge");
        require(Compound("[Na+].[Cl-]").charge() == 0, "Disconnected ions");
        require(Compound("c1ccccc1").formula().to_string() == "C6H6", "Aromatic formula");
        require(!Compound("F[C@H](Cl)Br").equals(Compound("F[C@@H](Cl)Br")), "Stereochemistry");
        const auto mapped = Compound("[OH:30][CH2:20][NH2:10]");
        require(mapped.smiles() == "NCO", "Map-independent canonicalization");
        require(mapped.get_atom_index_from_map(10) == 0, "Preserved map");
        require(mapped.get_atom_index_from_map(999) == -1, "Missing map");
        require(mapped.atom_map_to_index().at(30) == 2, "Map-to-index mapping");
        require(mapped.atom_index_mapped_smiles().find(":1") != std::string::npos, "Index maps");
        require(mapped.copy().equals(mapped), "Copy structure");
        require(mapped.copy().get_atom_index_from_map(30) == -1, "Copy regenerates maps");
        const auto restored = Compound::from_binary(mapped.mol_binary(true));
        require(restored.mapped_smiles() == mapped.mapped_smiles(), "Mapped binary round trip");
        require(Compound::from_binary(ethanol.mol_binary()).equals(ethanol), "Binary round trip");
        require(Compound("[CH3:50]O", true).get_atom_index_from_map(50) == -1, "Overwrite maps");
        require(Compound("[CH3:1]O").atom_map_to_index().size() == 2, "Fill missing maps");
        std::map<int, int> changes;
        const auto remapped = Compound::assign_atom_maps(mapped.mol_binary(true), true, changes);
        require(changes.size() == 3 && !remapped.empty(), "Map reassignment");
        rejects([&] { Compound::assign_atom_maps(mapped.mol_binary(true), true, changes); });
        rejects([] { Compound("invalid_smiles"); });
        rejects([] { Compound("[CH3:1][OH:1]"); });
        rejects([] { Compound::from_binary("invalid binary"); });
        std::cout << "Compound tests passed\n";
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
