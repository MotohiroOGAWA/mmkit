API Reference
=============

The ``mmkit`` API provides lightweight utilities for representing molecular
formulas and adduct ions in a consistent and programmatic way.

This reference documents the two core modules of the package:

- ``Formula``: tools for handling elemental compositions, exact masses, and charge-aware formula objects
- ``Adduct``: tools for representing adduct ions, parsing adduct strings, and converting neutral molecules into observed ion forms

These classes are designed to support mass spectrometry workflows such as
formula manipulation, adduct handling, and m/z calculation.

Classes
=======

.. list-table::
   :header-rows: 1

   * - Name
     - Description
   * - :class:`~mmkit.Formula.Formula`
     - Represents a molecular formula with elemental composition, charge, and exact mass.
   * - :class:`~mmkit.Adduct.Adduct`
     - Represents an adduct ion and provides conversion from neutral molecules to ion forms.


.. toctree::
   :maxdepth: 2
   :hidden:

   api/formula
   api/adduct