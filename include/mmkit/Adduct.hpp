#pragma once

#include "mmkit/Formula.hpp"
#include <utility>

namespace mmkit {

class Adduct {
public:
    Adduct(const std::string& ion_type = "M", int n_molecules = 1,
           const std::map<std::string, int>& formula_counts = {}, int charge = 0);
    static Adduct parse(const std::string& text);
    std::string ion_type() const;
    int n_molecules() const;
    int charge() const;
    std::map<std::string, int> formula_counts() const;
    Formula formula_diff() const;
    double mass_shift() const;
    std::map<std::string, int> element_diff() const;
    int get_formula_count(const std::string& formula) const;
    std::string to_string() const;
    Adduct copy() const;
    Adduct with_charge(int charge) const;
    bool equals(const Adduct& other) const;
    Adduct add(const Adduct& other, bool prefer_ion_type = false,
               bool prefer_n_molecules = false, bool prefer_charge = false) const;
    Adduct add_prefer_self(const Adduct& other) const;
    Formula apply_to_formula(const Formula& neutral) const;
    double apply_to_mass(double mass) const;
    double apply_to_mz(double mass) const;
    std::vector<Adduct> split(bool split_each = false) const;
    std::vector<bool> match_references(const std::vector<Adduct>& references) const;
    Adduct residual(const std::vector<Adduct>& references) const;
private:
    std::string ion_type_;
    int n_molecules_;
    std::vector<std::pair<Formula, int>> counts_;
    int charge_;
    std::pair<std::vector<bool>, Adduct> split_references(const std::vector<Adduct>& refs) const;
};
}
