#pragma once

#include "mmkit/Formula.hpp"
#include <memory>

namespace mmkit {

class Compound {
public:
    explicit Compound(const std::string& smiles, bool overwrite_atom_map = false);
    static Compound from_smiles(const std::string& smiles, bool overwrite_atom_map = false);
    // RDKit binary pickle bridge for Python. Not a portable storage format.
    static Compound from_binary(const std::string& binary, bool overwrite_atom_map = false);
    std::string mol_binary(bool mapped = false) const;
    static std::string assign_atom_maps(const std::string& binary, bool overwrite,
                                       std::map<int, int>& old_to_new);
    std::string smiles() const;
    std::string mapped_smiles() const;
    std::string atom_index_mapped_smiles() const;
    Formula formula() const;
    std::map<int, int> atom_map_to_index() const;
    int charge() const;
    double exact_mass() const;
    Compound copy() const;
    // -1 means no matching atom.
    int get_atom_index_from_map(int map_num) const;
    bool equals(const Compound& other) const;
private:
    struct Data;
    std::shared_ptr<const Data> data_;
    explicit Compound(std::shared_ptr<const Data> data);
};

}
