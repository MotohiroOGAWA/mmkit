#include "mmkit/Formula.hpp"

#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <regex>
#include <stdexcept>

namespace mmkit {
namespace {
int checked(long long value) {
    if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
        throw std::overflow_error("Count or charge exceeds 32-bit range");
    return static_cast<int>(value);
}
int integer(const std::string& text, bool negative) {
    try {
        const auto value = text.empty() ? 1LL : std::stoll(text);
        return checked(negative ? -value : value);
    } catch (const std::out_of_range&) {
        throw std::overflow_error("Count or charge exceeds 32-bit range");
    }
}
std::string charge_text(int charge) {
    if (!charge) return "";
    return std::string(charge > 0 ? "+" : "-") +
           (std::abs(static_cast<long long>(charge)) == 1 ? "" :
            std::to_string(std::abs(static_cast<long long>(charge))));
}
std::unique_ptr<RDKit::ROMol> read_smiles(const std::string& text) {
    try {
        auto mol = std::unique_ptr<RDKit::ROMol>(RDKit::SmilesToMol(text));
        if (!mol) throw std::invalid_argument("Invalid SMILES: " + text);
        return mol;
    } catch (const std::exception& e) {
        throw std::invalid_argument(std::string("Invalid SMILES: ") + e.what());
    }
}
Formula mol_formula(const RDKit::ROMol& mol) {
    std::unique_ptr<RDKit::ROMol> with_h(RDKit::MolOps::addHs(mol));
    std::map<std::string, int> counts;
    int charge = 0;
    for (const auto atom : with_h->atoms()) {
        ++counts[atom->getSymbol()];
        charge = checked(static_cast<long long>(charge) + atom->getFormalCharge());
    }
    return Formula(counts, charge);
}

}

Formula::Formula(const std::map<std::string, int>& elements, int charge, const std::string& raw)
    : charge_(charge), raw_(raw) {
    for (const auto& [symbol, count] : elements) {
        if (!count) continue;
        unsigned int number = 0;
        try { number = RDKit::PeriodicTable::getTable()->getAtomicNumber(symbol); }
        catch (const std::exception&) { throw std::invalid_argument("Unknown element: " + symbol); }
        if (!number) throw std::invalid_argument("Unknown element: " + symbol);
        elements_[symbol] = count;
    }
}
Formula Formula::parse(const std::string& text, bool store_raw) {
    std::string body = text;
    int charge = 0;
    std::smatch match;
    static const std::regex suffix("([+-])([0-9]*)$");
    if (std::regex_search(body, match, suffix)) {
        charge = integer(match[2], match[1] == "-");
        body = body.substr(0, match.position());
    }
    std::map<std::string, int> elements;
    static const std::regex term("([+-]?)([A-Z][a-z]?)([0-9]*)");
    auto begin = body.cbegin();
    while (begin != body.cend()) {
        std::match_results<std::string::const_iterator> part;
        if (!std::regex_search(begin, body.cend(), part, term, std::regex_constants::match_continuous))
            throw std::invalid_argument("Invalid formula: " + text);
        int value = integer(part[3], part[1] == "-");
        auto& count = elements[part[2]];
        count = checked(static_cast<long long>(count) + value);
        begin = part[0].second;
    }
    return Formula(elements, charge, store_raw ? text : "");
}
Formula Formula::from_smiles(const std::string& text) { return mol_formula(*read_smiles(text)); }
Formula Formula::from_binary(const std::string& binary) {
    RDKit::ROMol mol;
    RDKit::MolPickler::molFromPickle(binary, mol);
    return mol_formula(mol);
}
Formula Formula::empty() { return Formula(); }
std::map<std::string, int> Formula::elements() const { return elements_; }
std::vector<std::string> Formula::element_order() const {
    std::vector<std::string> keys;
    if (elements_.count("C")) {
        keys.push_back("C");
        if (elements_.count("H")) keys.push_back("H");
    }
    for (const auto& [symbol, count] : elements_)
        if (std::find(keys.begin(), keys.end(), symbol) == keys.end()) keys.push_back(symbol);
    return keys;
}
int Formula::charge() const { return charge_; }
std::string Formula::stored_raw_formula() const { return raw_; }
std::string Formula::raw_formula() const { return raw_.empty() ? to_string() : raw_; }
double Formula::exact_mass() const {
    double mass = 0;
    for (const auto& [symbol, count] : elements_)
        mass += RDKit::PeriodicTable::getTable()->getMostCommonIsotopeMass(symbol) * count;
    return mass;
}
bool Formula::is_nonnegative() const {
    return std::all_of(elements_.begin(), elements_.end(), [](const auto& e) { return e.second >= 0; });
}
std::string Formula::to_string(bool no_charge) const {
    std::string result;
    for (const auto& symbol : element_order()) {
        long long count = elements_.at(symbol);
        if (count < 0) result += "-";
        result += symbol;
        if (std::abs(count) != 1) result += std::to_string(std::abs(count));
    }
    return result + (no_charge ? "" : charge_text(charge_));
}
Formula Formula::plain() const { return with_charge(0); }
Formula Formula::normalized() const { return with_raw_formula(""); }
Formula Formula::normalized_plain() const { return Formula(elements_); }
Formula Formula::with_charge(int charge) const { return Formula(elements_, charge, raw_); }
Formula Formula::with_raw_formula(const std::string& raw) const { return Formula(elements_, charge_, raw); }
Formula Formula::copy() const { return *this; }
bool Formula::equals(const Formula& other, bool ignore_raw) const {
    return elements_ == other.elements_ && charge_ == other.charge_ && (ignore_raw || raw_ == other.raw_);
}
Formula Formula::add(const Formula& other) const {
    auto counts = elements_;
    for (const auto& [symbol, count] : other.elements_)
        counts[symbol] = checked(static_cast<long long>(counts[symbol]) + count);
    return Formula(counts, checked(static_cast<long long>(charge_) + other.charge_));
}
Formula Formula::subtract(const Formula& other) const {
    auto counts = elements_;
    for (const auto& [symbol, count] : other.elements_)
        counts[symbol] = checked(static_cast<long long>(counts[symbol]) - count);
    return Formula(counts, checked(static_cast<long long>(charge_) - other.charge_));
}
Formula Formula::multiply(int factor) const {
    auto counts = elements_;
    for (auto& [symbol, count] : counts) count = checked(static_cast<long long>(count) * factor);
    return Formula(counts, checked(static_cast<long long>(charge_) * factor), raw_);
}
bool Formula::contains(const Formula& other) const {
    if (!is_nonnegative()) throw std::invalid_argument("Containment requires a non-negative formula");
    for (const auto& [symbol, count] : other.elements_) {
        auto it = elements_.find(symbol);
        if ((it == elements_.end() ? 0 : it->second) < count) return false;
    }
    return true;
}

}
