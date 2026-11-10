#pragma once

#include <map>
#include <string>
#include <vector>

namespace mmkit {

// Immutable values. Counts may be negative; masses exclude electron correction.
class Formula {
public:
    Formula(const std::map<std::string, int>& elements = {}, int charge = 0,
            const std::string& raw_formula = "");
    static Formula parse(const std::string& text, bool store_raw = false);
    static Formula from_smiles(const std::string& smiles);
    static Formula from_binary(const std::string& binary);
    static Formula empty();
    std::map<std::string, int> elements() const;
    std::vector<std::string> element_order() const;
    int charge() const;
    std::string raw_formula() const;
    std::string stored_raw_formula() const;
    double exact_mass() const;
    bool is_nonnegative() const;
    std::string to_string(bool no_charge = false) const;
    Formula plain() const;
    Formula normalized() const;
    Formula normalized_plain() const;
    Formula with_charge(int charge) const;
    Formula with_raw_formula(const std::string& raw) const;
    Formula copy() const;
    bool equals(const Formula& other, bool ignore_raw = true) const;
    Formula add(const Formula& other) const;
    Formula subtract(const Formula& other) const;
    Formula multiply(int factor) const;
    bool contains(const Formula& other) const;
private:
    std::map<std::string, int> elements_;
    int charge_;
    std::string raw_;
};

}
