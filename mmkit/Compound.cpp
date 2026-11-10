#include "mmkit/Compound.hpp"

#include <GraphMol/Atom.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <regex>
#include <set>
#include <stdexcept>

namespace mmkit {
namespace {
int checked(long long value) {
    if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
        throw std::overflow_error("Count or charge exceeds 32-bit range");
    return static_cast<int>(value);
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
void assign_maps(RDKit::ROMol& mol, bool overwrite, std::map<int, int>& mapping) {
    std::set<int> used;
    for (auto atom : mol.atoms())
        if (!overwrite && atom->getAtomMapNum() > 0) used.insert(atom->getAtomMapNum());
    int next = 1;
    for (auto atom : mol.atoms()) {
        if (!overwrite && atom->getAtomMapNum()) continue;
        while (used.count(next)) next = checked(static_cast<long long>(next) + 1);
        int old = atom->getAtomMapNum();
        atom->setAtomMapNum(next);
        if (old > 0) mapping[old] = next;
        used.insert(next);
        next = checked(static_cast<long long>(next) + 1);
    }
}
}

struct Compound::Data {
    RDKit::ROMol mol;
    explicit Data(const RDKit::ROMol& source, bool overwrite) : mol(source) {
        std::vector<int> original;
        for (auto atom : mol.atoms()) {
            original.push_back(atom->getAtomMapNum());
            atom->setAtomMapNum(0);
        }
        const auto canonical = RDKit::MolToSmiles(mol);
        std::vector<unsigned int> order;
        mol.getProp("_smilesAtomOutputOrder", order);
        auto parsed = read_smiles(canonical);
        if (parsed->getNumAtoms() != order.size())
            throw std::invalid_argument("Atom count mismatch after canonicalization");
        if (!overwrite)
            for (unsigned int i = 0; i < order.size(); ++i)
                parsed->getAtomWithIdx(i)->setAtomMapNum(original.at(order[i]));
        std::map<int, int> ignored;
        assign_maps(*parsed, overwrite, ignored);
        std::set<int> seen;
        for (auto atom : parsed->atoms())
            if (!seen.insert(atom->getAtomMapNum()).second)
                throw std::invalid_argument("Duplicate atom map numbers");
        mol = std::move(*parsed);
    }
};
Compound::Compound(const std::string& text, bool overwrite)
    : data_(std::make_shared<Data>(*read_smiles(text), overwrite)) {}
Compound::Compound(std::shared_ptr<const Data> data) : data_(std::move(data)) {}
Compound Compound::from_smiles(const std::string& text, bool overwrite) { return Compound(text, overwrite); }
Compound Compound::from_binary(const std::string& binary, bool overwrite) {
    RDKit::ROMol mol;
    RDKit::MolPickler::molFromPickle(binary, mol);
    return Compound(std::make_shared<Data>(mol, overwrite));
}
std::string Compound::mol_binary(bool mapped) const {
    RDKit::ROMol mol(data_->mol);
    if (!mapped) for (auto atom : mol.atoms()) atom->setAtomMapNum(0);
    std::string binary;
    RDKit::MolPickler::pickleMol(mol, binary);
    return binary;
}
std::string Compound::assign_atom_maps(const std::string& binary, bool overwrite, std::map<int, int>& mapping) {
    if (!mapping.empty()) throw std::invalid_argument("Mapping must be empty");
    RDKit::ROMol mol;
    RDKit::MolPickler::molFromPickle(binary, mol);
    assign_maps(mol, overwrite, mapping);
    std::string result;
    RDKit::MolPickler::pickleMol(mol, result);
    return result;
}
std::string Compound::smiles() const {
    RDKit::ROMol mol(data_->mol);
    for (auto atom : mol.atoms()) atom->setAtomMapNum(0);
    return RDKit::MolToSmiles(mol);
}
std::string Compound::mapped_smiles() const {
    RDKit::ROMol mol(data_->mol);
    return RDKit::MolToSmiles(mol);
}
std::string Compound::atom_index_mapped_smiles() const {
    auto mol = read_smiles(smiles());
    for (auto atom : mol->atoms()) atom->setAtomMapNum(atom->getIdx());
    return RDKit::MolToSmiles(*mol);
}
Formula Compound::formula() const { return Formula::from_binary(mol_binary()); }
std::map<int, int> Compound::atom_map_to_index() const {
    RDKit::ROMol mol(data_->mol);
    for (auto atom : mol.atoms()) atom->setAtomMapNum(0);
    RDKit::MolToSmiles(mol);
    std::vector<unsigned int> order;
    mol.getProp("_smilesAtomOutputOrder", order);
    std::map<int, int> result;
    for (unsigned int i = 0; i < order.size(); ++i)
        result[data_->mol.getAtomWithIdx(order[i])->getAtomMapNum()] = i;
    return result;
}
int Compound::charge() const { return formula().charge(); }
double Compound::exact_mass() const { return formula().exact_mass(); }
Compound Compound::copy() const { return Compound::from_binary(mol_binary(true), true); }
int Compound::get_atom_index_from_map(int number) const {
    for (auto atom : data_->mol.atoms())
        if (atom->getAtomMapNum() == number) return atom->getIdx();
    return -1;
}
bool Compound::equals(const Compound& other) const { return smiles() == other.smiles(); }

}
