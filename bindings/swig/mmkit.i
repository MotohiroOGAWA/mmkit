%module mmkit
%{
#include "mmkit/mmkit.hpp"
%}
%include <std_string.i>
%include <std_map.i>
%include <std_vector.i>
%include <exception.i>
%exception {
  try { $action }
  catch (const std::invalid_argument& e) { SWIG_exception(SWIG_ValueError, e.what()); }
  catch (const std::out_of_range& e) { SWIG_exception(SWIG_IndexError, e.what()); }
  catch (const std::exception& e) { SWIG_exception(SWIG_RuntimeError, e.what()); }
  catch (...) { SWIG_exception(SWIG_RuntimeError, "Unknown native error"); }
}
%ignore mmkit::Formula::from_binary;
%ignore mmkit::Compound::from_binary;
%ignore mmkit::Compound::mol_binary;
%ignore mmkit::Compound::assign_atom_maps;
#ifdef SWIGJAVA
%typemap(javainterfaces) mmkit::Formula "java.lang.AutoCloseable";
%typemap(javainterfaces) mmkit::Compound "java.lang.AutoCloseable";
%typemap(javainterfaces) mmkit::Adduct "java.lang.AutoCloseable";
%typemap(javacode) mmkit::Formula %{ public void close() { delete(); } %}
%typemap(javacode) mmkit::Compound %{ public void close() { delete(); } %}
%typemap(javacode) mmkit::Adduct %{ public void close() { delete(); } %}
#endif
namespace std {
  %template(StringIntMap) map<string, int>;
  %template(IntIntMap) map<int, int>;
  %template(StringVector) vector<string>;
  %template(BoolVector) vector<bool>;
}
%include "mmkit/Formula.hpp"
%include "mmkit/Compound.hpp"
%include "mmkit/Adduct.hpp"
namespace std {
  %template(AdductVector) vector<mmkit::Adduct>;
}
