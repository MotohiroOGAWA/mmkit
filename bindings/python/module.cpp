#include "mmkit/mmkit.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <sstream>
namespace py = pybind11;
using namespace mmkit;

namespace {
py::object chem() { return py::module_::import("rdkit.Chem"); }
std::string mol_bytes(py::object mol) {
    if (!py::isinstance(mol, chem().attr("Mol"))) throw py::type_error("mol must be an RDKit Mol object");
    return mol.attr("ToBinary")().cast<std::string>();
}
py::object as_mol(const std::string& binary) { return chem().attr("Mol")(py::bytes(binary)); }
std::string formula_key(py::handle key) {
    if (py::isinstance<Formula>(key)) return key.cast<Formula>().plain().to_string();
    if (py::isinstance<py::str>(key)) return key.cast<std::string>();
    throw py::type_error("formula must be str or Formula");
}
std::map<std::string, int> formula_counts(py::dict counts) {
    std::map<std::string, int> result;
    for (auto item : counts) result[formula_key(item.first)] = py::cast<int>(item.second);
    return result;
}
py::dict adduct_formulas(const Adduct& value) {
    py::dict result;
    for (const auto& [key, count] : value.formula_counts()) result[py::cast(Formula::parse(key, true).plain())] = count;
    return result;
}
py::tuple as_tuple(const std::vector<Adduct>& values) {
    py::tuple result(values.size());
    for (std::size_t i = 0; i < values.size(); ++i) result[i] = py::cast(values[i]);
    return result;
}
}
PYBIND11_MODULE(_native, m) {
    m.doc() = "RDKit-backed mmkit C++ core";
    py::class_<Formula>(m, "Formula")
        .def(py::init<const std::map<std::string,int>&,int,const std::string&>(), py::arg("elements"), py::arg("charge")=0, py::arg("raw_formula")="")
        .def_static("parse", &Formula::parse, py::arg("formula_str"), py::arg("store_raw")=false)
        .def_static("from_dict", [](const std::map<std::string,int>& e,int c,const std::string& r) { return Formula(e,c,r); }, py::arg("elements"), py::kw_only(), py::arg("charge")=0, py::arg("raw_formula")="")
        .def_static("from_smiles", &Formula::from_smiles)
        .def_static("from_mol", [](py::object mol) { return Compound::from_binary(mol_bytes(mol)).formula(); })
        .def_static("empty", &Formula::empty)
        .def_property_readonly("elements", [](const Formula& f) {
            py::object result = py::module_::import("collections").attr("OrderedDict")();
            const auto counts = f.elements();
            for (const auto& symbol : f.element_order()) result[py::str(symbol)] = counts.at(symbol);
            return result;
        })
        .def_property_readonly("charge", &Formula::charge)
        .def_property_readonly("raw_formula", &Formula::raw_formula)
        .def_property_readonly("_raw_formula", &Formula::stored_raw_formula)
        .def_property_readonly("exact_mass", &Formula::exact_mass)
        .def_property_readonly("is_nonnegative", &Formula::is_nonnegative)
        .def_property_readonly("value", [](const Formula& f) { return f.to_string(); })
        .def_property_readonly("plain_value", [](const Formula& f) { return f.to_string(true); })
        .def_property_readonly("plain", &Formula::plain)
        .def_property_readonly("normalized", &Formula::normalized)
        .def_property_readonly("normalized_plain", &Formula::normalized_plain)
        .def("with_charge", &Formula::with_charge, py::arg("charge"))
        .def("with_raw_formula", &Formula::with_raw_formula, py::arg("raw_formula"))
        .def("copy", &Formula::copy)
        .def("to_dict", &Formula::elements)
        .def("to_string", &Formula::to_string, py::arg("no_charge")=false)
        .def("eq", &Formula::equals, py::arg("other"), py::arg("ignore_raw")=true)
        .def("__str__", [](const Formula& f) { return f.to_string(); })
        .def("__repr__", [](const Formula& f) { return "Formula(" + f.to_string() + ")"; })
        .def("__eq__", [](const Formula& f, py::object o) { return py::isinstance<Formula>(o) && f.equals(o.cast<Formula>()); }, py::is_operator())
        .def("__hash__", [](const Formula& f) { return py::hash(py::str(f.to_string())); })
        .def("__add__", &Formula::add, py::is_operator())
        .def("__sub__", &Formula::subtract, py::is_operator())
        .def("__mul__", &Formula::multiply, py::is_operator())
        .def("__rmul__", &Formula::multiply, py::is_operator())
        .def("__contains__", [](const Formula& f, py::object o) { return f.contains(Formula::parse(formula_key(o))); });

    py::class_<Compound>(m, "Compound")
        .def(py::init([](py::object mol,bool overwrite) { return Compound::from_binary(mol_bytes(mol),overwrite); }), py::arg("mol"), py::arg("overwrite_atom_map")=false)
        .def_static("from_smiles", &Compound::from_smiles, py::arg("smiles"), py::arg("overwrite_atom_map")=false)
        .def_property_readonly("mol", [](const Compound& c) { return as_mol(c.mol_binary()); })
        .def_property_readonly("mapped_mol", [](const Compound& c) { return as_mol(c.mol_binary(true)); })
        .def_property_readonly("smiles", &Compound::smiles)
        .def_property_readonly("mapped_smiles", &Compound::mapped_smiles)
        .def_property_readonly("atom_index_mapped_smiles", &Compound::atom_index_mapped_smiles)
        .def_property_readonly("formula", &Compound::formula)
        .def_property_readonly("charge", &Compound::charge)
        .def_property_readonly("exact_mass", &Compound::exact_mass)
        .def_property_readonly("atom_map_to_index", [](const Compound& c) { return py::module_::import("bidict").attr("bidict")(py::cast(c.atom_map_to_index())); })
        .def("copy", &Compound::copy)
        .def("get_atom_index_from_map", [](const Compound& c,int n) -> py::object {
            int i = c.get_atom_index_from_map(n);
            if (i < 0) return py::none();
            return py::int_(i);
        }, py::arg("map_num"))
        .def_static("_smiles_atom_output_order", [](py::object mol) {
            std::string text = mol.attr("GetProp")("_smilesAtomOutputOrder").cast<std::string>();
            for (auto& c : text) if (c == '[' || c == ']' || c == ',') c = ' ';
            std::istringstream input(text);
            std::vector<int> order;
            int index;
            while (input >> index) order.push_back(index);
            return order;
        })
        .def_static("_assign_atom_map", [](py::object mol, bool overwrite, py::object mapping) {
            if (!mapping.is_none() && !py::isinstance<py::dict>(mapping)) throw py::type_error("atom_map_dict must be a dict");
            if (!mapping.is_none() && py::len(mapping)) throw py::value_error("atom_map_dict must be empty");
            std::map<int,int> changes;
            auto result = Compound::assign_atom_maps(mol_bytes(mol),overwrite,changes);
            if (!mapping.is_none()) for (const auto& [a,b] : changes) mapping[py::int_(a)] = b;
            return as_mol(result);
        }, py::arg("mol"), py::kw_only(), py::arg("overwrite")=false, py::arg("atom_map_dict")=py::none())
        .def("__str__", &Compound::smiles)
        .def("__repr__", [](const Compound& c) { return "Compound(smiles=" + c.smiles() + ")"; })
        .def("__eq__", [](const Compound& c,py::object o) { return py::isinstance<Compound>(o) && c.equals(o.cast<Compound>()); },py::is_operator())
        .def("__hash__", [](const Compound& c) { return py::hash(py::str(c.smiles())); });

    py::class_<Adduct>(m, "Adduct")
        .def(py::init([](const std::string& type,int n,py::dict counts,int charge) { return Adduct(type,n,formula_counts(counts),charge); }), py::arg("ion_type")="M",py::arg("n_molecules")=1,py::arg("formula_counts")=py::dict(),py::arg("charge")=0)
        .def_static("parse", &Adduct::parse, py::arg("adduct_str"))
        .def_static("from_dict", [](py::dict counts,const std::string& type,int n,int charge) { return Adduct(type,n,formula_counts(counts),charge); }, py::arg("formula_counts"),py::kw_only(),py::arg("ion_type")="M",py::arg("n_molecules")=1,py::arg("charge")=0)
        .def_static("from_adducts", [](const std::string& type,int n,py::tuple in,py::tuple out,int charge) {
            auto result = Adduct(type,n,{},charge);
            for (auto item : in) result = result.add_prefer_self(Adduct(type,n,{{formula_key(item),1}},charge));
            for (auto item : out) result = result.add_prefer_self(Adduct(type,n,{{formula_key(item),-1}},charge));
            return result;
        }, py::kw_only(),py::arg("ion_type")="M",py::arg("n_molecules")=1,py::arg("adducts_in")=py::tuple(),py::arg("adducts_out")=py::tuple(),py::arg("charge")=0)
        .def_property_readonly("ion_type", &Adduct::ion_type)
        .def_property_readonly("n_molecules", &Adduct::n_molecules)
        .def_property_readonly("charge", &Adduct::charge)
        .def_property_readonly("formula_counts", &adduct_formulas)
        .def_property_readonly("adduct_formulas", &adduct_formulas)
        .def_property_readonly("formula_diff", &Adduct::formula_diff)
        .def_property_readonly("mass_shift", &Adduct::mass_shift)
        .def_property_readonly("element_diff", &Adduct::element_diff)
        .def("copy", &Adduct::copy)
        .def("with_charge", &Adduct::with_charge,py::arg("charge"))
        .def("get_formula_count", [](const Adduct& a,py::object key) { return a.get_formula_count(formula_key(key)); },py::arg("formula"))
        .def("to_dict", &Adduct::formula_counts)
        .def("add", &Adduct::add,py::arg("other"),py::kw_only(),py::arg("prefer_ion_type")=false,py::arg("prefer_n_molecules")=false,py::arg("prefer_charge")=false)
        .def("add_prefer_self", &Adduct::add_prefer_self,py::arg("other"))
        .def("apply_to_formula", &Adduct::apply_to_formula,py::arg("neutral_formula"))
        .def("apply_to_mass", &Adduct::apply_to_mass,py::arg("neutral_mass"))
        .def("apply_to_mz", &Adduct::apply_to_mz,py::arg("neutral_mass"))
        .def("split", [](const Adduct& a,bool each) { return as_tuple(a.split(each)); },py::arg("split_each")=false)
        .def_static("split_by_reference_adducts", [](const Adduct& a,const std::vector<Adduct>& refs) {
            auto flags = a.match_references(refs);
            py::tuple matches(flags.size());
            for (std::size_t i=0;i<flags.size();++i) matches[i]=py::bool_(flags[i]);
            return py::make_tuple(matches,a.residual(refs));
        },py::arg("adduct"),py::arg("reference_adducts"))
        .def("__str__", &Adduct::to_string)
        .def("__repr__", [](const Adduct& a) { return "Adduct(" + a.to_string() + ")"; })
        .def("__eq__", [](const Adduct& a,py::object o) { return py::isinstance<Adduct>(o) && a.equals(o.cast<Adduct>()); },py::is_operator())
        .def("__hash__", [](const Adduct& a) { return py::hash(py::str(a.to_string())); });
}
