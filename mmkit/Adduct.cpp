#include "mmkit/Adduct.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <regex>
#include <stdexcept>

namespace mmkit {
namespace {
int checked(long long value) {
    if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
        throw std::overflow_error("Count or charge exceeds 32-bit range");
    return static_cast<int>(value);
}
int integer(const std::string& text, int fallback = 1) {
    if (text.empty()) return fallback;
    return checked(std::stoll(text));
}
std::string charge_text(int charge) {
    if (!charge) return "";
    return std::string(charge > 0 ? "+" : "-") +
           (std::abs(static_cast<long long>(charge)) == 1 ? "" :
            std::to_string(std::abs(static_cast<long long>(charge))));
}
}

Adduct::Adduct(const std::string& type, int n, const std::map<std::string, int>& counts, int charge)
    : ion_type_(type), n_molecules_(n), charge_(charge) {
    if (type != "M" && type != "F") throw std::invalid_argument("ion_type must be M or F");
    if (n < 1) throw std::invalid_argument("n_molecules must be at least 1");
    for (const auto& [text, count] : counts) {
        if (!count) continue;
        auto formula = Formula::parse(text, true).plain();
        if (!formula.is_nonnegative() || formula.elements().empty())
            throw std::invalid_argument("Adduct terms must be nonempty, non-negative formulas");
        auto it = std::find_if(counts_.begin(), counts_.end(), [&](const auto& e) { return e.first.equals(formula); });
        if (it == counts_.end()) counts_.emplace_back(formula, count);
        else it->second = checked(static_cast<long long>(it->second) + count);
    }
    counts_.erase(std::remove_if(counts_.begin(), counts_.end(), [](const auto& e) { return !e.second; }), counts_.end());
}
Adduct Adduct::parse(const std::string& text) {
    static const std::regex format(R"(^\[([0-9]*)(M|F)([^\]]*)\]([+-][0-9]*|[0-9]+[+-])?$)");
    std::smatch parts;
    if (!std::regex_match(text, parts, format)) throw std::invalid_argument("Invalid adduct: " + text);
    std::string suffix = parts[4];
    int charge = 0;
    if (!suffix.empty()) {
        bool first = suffix.front() == '+' || suffix.front() == '-';
        char sign = first ? suffix.front() : suffix.back();
        charge = checked(static_cast<long long>(integer(first ? suffix.substr(1) : suffix.substr(0, suffix.size()-1))) * (sign == '-' ? -1 : 1));
    }
    std::map<std::string, int> counts;
    std::string remainder = parts[3];
    auto begin = remainder.cbegin();
    static const std::regex term("([+-])([0-9]*)([A-Z][A-Za-z0-9]*)");
    while (begin != remainder.cend()) {
        std::match_results<std::string::const_iterator> item;
        if (!std::regex_search(begin, remainder.cend(), item, term, std::regex_constants::match_continuous))
            throw std::invalid_argument("Invalid adduct term: " + text);
        auto& count = counts[item[3]];
        count = checked(static_cast<long long>(count) + static_cast<long long>(integer(item[2])) * (item[1] == "-" ? -1 : 1));
        begin = item[0].second;
    }
    return Adduct(parts[2], integer(parts[1]), counts, charge);
}
std::string Adduct::ion_type() const { return ion_type_; }
int Adduct::n_molecules() const { return n_molecules_; }
int Adduct::charge() const { return charge_; }
std::map<std::string, int> Adduct::formula_counts() const {
    std::map<std::string, int> result;
    for (const auto& [formula, count] : counts_) {
        auto name = formula.raw_formula();
        name.erase(std::remove(name.begin(), name.end(), '+'), name.end());
        name.erase(std::remove(name.begin(), name.end(), '-'), name.end());
        result[name] = count;
    }
    return result;
}
Formula Adduct::formula_diff() const {
    auto result = Formula::empty();
    for (const auto& [formula, count] : counts_) result = result.add(formula.multiply(count));
    return result;
}
double Adduct::mass_shift() const {
    double result = 0;
    for (const auto& [formula, count] : counts_) result += formula.exact_mass() * count;
    return result;
}
std::map<std::string, int> Adduct::element_diff() const {
    // Retain zero totals for compatibility with the Python API.
    std::map<std::string, int> result;
    for (const auto& [formula, count] : counts_)
        for (const auto& [symbol, n] : formula.elements())
            result[symbol] = checked(static_cast<long long>(result[symbol]) + static_cast<long long>(n) * count);
    return result;
}
int Adduct::get_formula_count(const std::string& text) const {
    auto formula = Formula::parse(text).plain();
    for (const auto& item : counts_) if (item.first.equals(formula)) return item.second;
    return 0;
}
std::string Adduct::to_string() const {
    auto ordered = counts_;
    std::sort(ordered.begin(), ordered.end(), [](const auto& a, const auto& b) {
        double am = a.first.exact_mass(), bm = b.first.exact_mass();
        return am == bm ? a.first.raw_formula() < b.first.raw_formula() : am > bm;
    });
    std::string result = "[" + (n_molecules_ > 1 ? std::to_string(n_molecules_) : "") + ion_type_;
    for (const auto& [formula, count] : ordered) {
        long long n = count;
        result += n > 0 ? "+" : "-";
        if (std::abs(n) != 1) result += std::to_string(std::abs(n));
        auto text = formula.raw_formula();
        text.erase(std::remove(text.begin(), text.end(), '+'), text.end());
        text.erase(std::remove(text.begin(), text.end(), '-'), text.end());
        result += text;
    }
    return result + "]" + charge_text(charge_);
}
Adduct Adduct::copy() const { return *this; }
Adduct Adduct::with_charge(int charge) const { auto result = *this; result.charge_ = charge; return result; }
bool Adduct::equals(const Adduct& other) const { return to_string() == other.to_string(); }
Adduct Adduct::add(const Adduct& other, bool type, bool n, bool charge) const {
    if ((!type && ion_type_ != other.ion_type_) || (!n && n_molecules_ != other.n_molecules_) ||
        (!charge && charge_ != other.charge_)) throw std::invalid_argument("Incompatible adduct metadata");
    auto counts = formula_counts();
    for (const auto& [text, count] : other.formula_counts()) counts[text] = checked(static_cast<long long>(counts[text]) + count);
    return Adduct(ion_type_, n_molecules_, counts, charge_);
}
Adduct Adduct::add_prefer_self(const Adduct& other) const { return add(other, true, true, true); }
Formula Adduct::apply_to_formula(const Formula& neutral) const {
    return neutral.multiply(n_molecules_).add(formula_diff()).with_charge(charge_);
}
double Adduct::apply_to_mass(double mass) const {
    if (!std::isfinite(mass)) throw std::invalid_argument("Mass must be finite");
    return mass * n_molecules_ + mass_shift();
}
double Adduct::apply_to_mz(double mass) const {
    if (!charge_) throw std::invalid_argument("Cannot calculate m/z for an uncharged adduct");
    return apply_to_mass(mass) / std::abs(static_cast<double>(charge_));
}
std::vector<Adduct> Adduct::split(bool each) const {
    std::vector<Adduct> result;
    for (int sign : {1, -1}) {
        std::map<std::string, int> group;
        for (const auto& [formula, count] : counts_) {
            if ((count > 0 ? 1 : -1) != sign) continue;
            if (each) {
                for (long long i = 0; i < std::abs(static_cast<long long>(count)); ++i)
                    result.emplace_back(ion_type_, n_molecules_, std::map<std::string, int>{{formula.raw_formula(), sign}}, 0);
            } else group[formula.raw_formula()] = count;
        }
        if (!each && !group.empty()) result.emplace_back(ion_type_, n_molecules_, group, 0);
    }
    return result;
}
std::pair<std::vector<bool>, Adduct> Adduct::split_references(const std::vector<Adduct>& refs) const {
    if (!charge_) throw std::invalid_argument("Adduct must have non-zero charge");
    auto remaining = *this;
    std::vector<bool> flags;
    for (const auto& ref : refs) {
        bool match = ref.ion_type_ == ion_type_ && ref.n_molecules_ == n_molecules_ &&
                     ((ref.charge_ > 0 && charge_ > 0) || (ref.charge_ < 0 && charge_ < 0));
        for (const auto& [formula, count] : ref.counts_)
            if (remaining.get_formula_count(formula.to_string(true)) != count) match = false;
        flags.push_back(match);
        if (match) {
            for (const auto& [formula, count] : ref.counts_)
                for (auto& item : remaining.counts_)
                    if (item.first.equals(formula)) item.second = checked(static_cast<long long>(item.second) - count);
        }
    }
    remaining.counts_.erase(std::remove_if(remaining.counts_.begin(), remaining.counts_.end(), [](const auto& e) { return !e.second; }), remaining.counts_.end());
    remaining.n_molecules_ = 1;
    remaining.charge_ = 0;
    return {flags, remaining};
}
std::vector<bool> Adduct::match_references(const std::vector<Adduct>& refs) const { return split_references(refs).first; }
Adduct Adduct::residual(const std::vector<Adduct>& refs) const { return split_references(refs).second; }
}
