#include "mmkit/mmkit.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

using mmkit::Adduct;
using mmkit::Formula;
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void near(double actual, double expected) {
    require(std::abs(actual - expected) < 1e-6, "Mass mismatch");
}
template<class Exception = std::invalid_argument, class Action> void rejects(Action action) {
    try { action(); } catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception");
}
int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"construction and parsing", [] {
            require(Adduct().to_string() == "[M]", "Default");
            for (const auto& text : {"[M+H]+", "[M-H]-", "[2M+Na]+", "[F+2H]+2", "[M-5H]-"})
                require(Adduct::parse(text).to_string() == text, "Round trip");
            require(Adduct::parse("[M+2H]2+").equals(Adduct::parse("[M+2H]+2")), "Charge suffix");
            const auto a = Adduct::parse("[M+HCOOH-H]-");
            require(a.ion_type() == "M" && a.n_molecules() == 1 && a.charge() == -1, "Metadata");
            require(a.get_formula_count("CH2O2") == 1 && a.get_formula_count("H+") == -1, "Formula lookup");
            require(a.get_formula_count("Na") == 0, "Absent formula");
            require(Adduct::parse("[M+H+H-H]+").to_string() == "[M+H]+", "Repeated terms");
            require(Adduct::parse("[M+H-H]").to_string() == "[M]", "Cancellation");
            require(Adduct("M", 1, {{"HOH", 1}, {"H2O", -1}}).to_string() == "[M]", "Equivalent formula cancellation");
        }},
        {"value semantics and composition", [] {
            const auto a = Adduct::parse("[M+H]+");
            require(a.copy().equals(a), "Copy");
            require(a.with_charge(-2).charge() == -2 && a.charge() == 1, "Immutable charge");
            auto counts = a.formula_counts(); counts["H"] = 99;
            require(a.get_formula_count("H") == 1, "Independent counts");
            require(a.formula_diff().equals(Formula::parse("H")), "Formula difference");
            near(a.mass_shift(), Formula::parse("H").exact_mass());
            const auto loss = Adduct::parse("[M+H-H2O]+");
            require(loss.element_diff().at("H") == -1 && loss.element_diff().at("O") == -1, "Element difference");
            require(a.add(Adduct::parse("[M+Na]+")).to_string() == "[M+Na+H]+", "Addition");
            require(a.add_prefer_self(Adduct::parse("[2F-H]-")).to_string() == "[M]+", "Prefer self");
            rejects([&] { a.add(Adduct("F")); });
            rejects([&] { a.add(Adduct::parse("[2M+H]+")); });
            rejects([&] { a.add(Adduct::parse("[M-H]-")); });
        }},
        {"Compound and Formula integration", [] {
            const auto compound = mmkit::Compound::from_smiles("CCO");
            const auto protonated = Adduct::parse("[M+H]+");
            require(protonated.apply_to_formula(compound.formula()).to_string() == "C2H7O+", "Protonated formula");
            near(protonated.apply_to_mz(compound.exact_mass()), 47.049689812);
            const auto dimer = Adduct::parse("[2M+2H]+2");
            require(dimer.apply_to_formula(compound.formula()).to_string() == "C4H14O2+2", "Dimer formula");
            near(dimer.apply_to_mz(compound.exact_mass()), protonated.apply_to_mz(compound.exact_mass()));
            const auto negative = Adduct::parse("[M-2H]-2");
            near(negative.apply_to_mz(100), (100 - 2 * Formula::parse("H").exact_mass()) / 2);
            near(Adduct().apply_to_mass(100), 100);
            rejects([] { Adduct().apply_to_mz(100); });
            rejects([] { Adduct().apply_to_mass(std::numeric_limits<double>::infinity()); });
            rejects([] { Adduct().apply_to_mass(std::numeric_limits<double>::quiet_NaN()); });
        }},
        {"splitting and references", [] {
            const auto a = Adduct::parse("[M+2H+Na-H2O]+");
            const auto grouped = a.split();
            require(grouped.size() == 2 && grouped[0].charge() == 0 && grouped[1].get_formula_count("H2O") == -1, "Grouped split");
            require(a.split(true).size() == 4, "Individual split");
            require(Adduct().split().empty(), "Empty split");
            const std::vector<Adduct> refs = {Adduct::parse("[M+H]+"), Adduct::parse("[M+2H]+"), Adduct::parse("[M+Na]+"), Adduct::parse("[M-H]-")};
            require(a.match_references(refs) == std::vector<bool>({false, true, true, false}), "Exact counts and polarity");
            require(a.residual(refs).to_string() == "[M-H2O]", "Residual loss");
            require(a.residual({}).with_charge(1).equals(a), "No references");
            const auto dimer = Adduct::parse("[2M+H-H2O]+");
            require(dimer.match_references({refs[0]}) == std::vector<bool>({false}), "Molecule count mismatch");
            require(dimer.residual({Adduct::parse("[2M+H]+")}).to_string() == "[M-H2O]", "Dimer residual");
            require(Adduct::parse("[M-H]-").match_references({refs[0], refs[3]}) == std::vector<bool>({false, true}), "Negative reference");
            require(Adduct::parse("[M+H]+").match_references({refs[0], refs[0]}) == std::vector<bool>({true, false}), "Consume once");
            rejects([] { Adduct().match_references({}); });
        }},
        {"invalid inputs and overflow", [] {
            for (const auto& text : {"", "M+H+", "abc", "[+H]", "[0M+H]+", "[M+H!]+", "[M+Xx]+", "[M+]", "[M+H]++"})
                rejects([&] { Adduct::parse(text); });
            rejects([] { Adduct("X"); });
            rejects([] { Adduct("M", 0); });
            rejects([] { Adduct("M", 1, {{"-H", 1}}); });
            rejects([] { Adduct("M", 1, {{"", 1}}); });
            rejects<std::overflow_error>([] { Adduct::parse("[M+2147483647H+H]+"); });
            rejects<std::overflow_error>([] { Adduct::parse("[M+2147483647H]+").add(Adduct::parse("[M+H]+")); });
            rejects<std::overflow_error>([] { Adduct::parse("[2M]").apply_to_formula(Formula({{"C", std::numeric_limits<int>::max()}})); });
        }}
    };
    int failures = 0;
    for (const auto& [name, test] : tests) {
        try { test(); std::cout << "PASS: " << name << '\n'; }
        catch (const std::exception& e) { ++failures; std::cerr << "FAIL: " << name << ": " << e.what() << '\n'; }
    }
    return failures ? 1 : 0;
}
