#include "mmkit/Formula.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using mmkit::Formula;
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
template<class Exception, class Action> void throws(Action action) {
    try { action(); } catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception was not thrown");
}
int main() {
    const std::vector<std::pair<std::string, std::function<void()>>> tests = {
        {"parse and Hill order", [] {
            for (const auto& item : std::vector<std::pair<std::string, std::string>>{
                {"C6O6H12", "C6H12O6"}, {"HOH", "H2O"},
                {"NaClH", "ClHNa"}, {"-H5-C-OH", "-C-H4-O"},
                {"C6H12-NaO6-", "C6H12-NaO6-"}, {"H0", ""},
                {"H-H", ""}, {"Na+2", "Na+2"}, {"O-2", "O-2"}})
                require(Formula::parse(item.first).to_string() == item.second, "Unexpected normalized formula");
            require(Formula::empty().equals(Formula::parse("")), "Empty formula");
        }},
        {"invalid input", [] {
            for (const auto& text : {"H2O!", " H2O", "h2o", "(OH)2", "2H", "H++", "Xx"})
                throws<std::invalid_argument>([&] { Formula::parse(text); });
            throws<std::invalid_argument>([] { Formula({{"Unknown", 1}}); });
            throws<std::invalid_argument>([] { Formula::from_smiles("invalid_smiles"); });
            throws<std::exception>([] { Formula::from_binary("invalid binary"); });
        }},
        {"raw and normalization", [] {
            const auto f = Formula::parse("OH2+", true);
            require(f.raw_formula() == "OH2+", "Raw formula lost");
            require(f.plain().to_string() == "H2O" && f.plain().raw_formula() == "OH2+", "Plain");
            require(f.normalized().stored_raw_formula().empty() && f.normalized().charge() == 1, "Normalized");
            require(f.normalized_plain().charge() == 0 && f.normalized_plain().stored_raw_formula().empty(), "Normalized plain");
            require(f.equals(Formula::parse("H2O+", true)), "Equality ignores raw");
            require(!f.equals(Formula::parse("H2O+", true), false), "Raw equality");
            require(f.copy().equals(f, false), "Copy");
            require(f.with_raw_formula("custom").raw_formula() == "custom", "Raw replacement");
            require(f.with_charge(-2).to_string() == "H2O-2" && f.charge() == 1, "Value immutability");
        }},
        {"arithmetic and containment", [] {
            const auto glucose = Formula::parse("C6H12O6");
            const auto water = Formula::parse("H2O");
            require(glucose.add(water).to_string() == "C6H14O7", "Addition");
            require(glucose.subtract(water).to_string() == "C6H10O5", "Subtraction");
            require(glucose.subtract(glucose).equals(Formula::empty()), "Cancellation");
            require(water.multiply(-2).to_string() == "-H4-O2", "Negative multiplication");
            require(water.multiply(0).equals(Formula::empty()), "Zero multiplication");
            require(Formula::parse("Na+").multiply(2).to_string() == "Na2+2", "Charge multiplication");
            require(glucose.contains(water) && glucose.contains(Formula::empty()), "Containment");
            require(!glucose.contains(Formula::parse("Na")), "Missing element");
            require(!glucose.contains(Formula::parse("H13")), "Insufficient count");
            throws<std::invalid_argument>([] { Formula::parse("-H").contains(Formula::parse("H")); });
            auto elements = water.elements(); elements["H"] = 99;
            require(water.to_string() == "H2O", "Elements must be independent");
        }},
        {"mass and molecules", [] {
            require(std::abs(Formula::parse("H2O").exact_mass() - 18.010565) < 0.000001, "Water mass");
            require(std::abs(Formula::parse("C6H12O6").exact_mass() - 180.063388) < 0.000001, "Glucose mass");
            require(Formula::empty().exact_mass() == 0, "Empty mass");
            require(Formula::parse("-H").exact_mass() < 0, "Negative mass");
            require(Formula::parse("Na+").exact_mass() == Formula::parse("Na").exact_mass(), "No electron correction");
            require(Formula::from_smiles("O").to_string() == "H2O", "Implicit hydrogens");
            require(Formula::from_smiles("[NH4+]").to_string() == "H4N+", "Formal charge");
            require(Formula::from_smiles("c1ccccc1").to_string() == "C6H6", "Aromatic hydrogens");
            require(Formula::from_smiles("[Na+].[Cl-]").to_string() == "ClNa", "Disconnected molecule");
        }},
        {"integer boundaries", [] {
            constexpr int hi = std::numeric_limits<int>::max();
            constexpr int lo = std::numeric_limits<int>::min();
            const Formula minimum({{"C", lo}}, lo);
            require(Formula::parse(minimum.to_string()).equals(minimum), "INT_MIN round trip");
            require(minimum.subtract(minimum).equals(Formula::empty()), "INT_MIN cancellation");
            require(Formula({{"C", hi}}, hi).subtract(Formula({{"C", hi}}, hi)).equals(Formula::empty()), "INT_MAX cancellation");
            throws<std::overflow_error>([&] { Formula({{"C", hi}}).add(Formula::parse("C")); });
            throws<std::overflow_error>([&] { Formula({}, hi).add(Formula({}, 1)); });
            throws<std::overflow_error>([&] { minimum.subtract(Formula::parse("C")); });
            throws<std::overflow_error>([&] { minimum.multiply(-1); });
            throws<std::overflow_error>([&] { Formula({}, lo).subtract(Formula({}, 1)); });
            for (const auto& text : {"C2147483648", "C+2147483648", "C999999999999999999999999999", "C2147483647C"})
                throws<std::overflow_error>([&] { Formula::parse(text); });
        }}
    };
    int failures = 0;
    for (const auto& [name, test] : tests) {
        try { test(); std::cout << "PASS: " << name << '\n'; }
        catch (const std::exception& e) { ++failures; std::cerr << "FAIL: " << name << ": " << e.what() << '\n'; }
    }
    return failures ? 1 : 0;
}
