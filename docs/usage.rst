Usage
=====

Formula parsing
---------------

You can parse a chemical formula string into a ``Formula`` object.

.. code-block:: python

   from mmkit.chem.Formula import Formula

   f = Formula.parse("C6H12O6")
   print(f.value)

Formula arithmetic
------------------

You can add, subtract, and multiply formulas.

.. code-block:: python

   water = Formula.parse("H2O")
   print((water * 2).value)
   print((Formula.parse("C6H12O6") + water).value)