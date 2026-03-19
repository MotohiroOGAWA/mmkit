{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

{% if docstring %}
{{ docstring }}
{% endif %}

{% if attributes %}
Properties
----------

.. autosummary::
   :nosignatures:

{% for name in attributes %}
   {{ objname }}.{{ name }}
{% endfor %}
{% endif %}

{% if methods %}
Methods
-------

.. autosummary::
   :nosignatures:

{% for name in methods %}
   {{ objname }}.{{ name }}
{% endfor %}
{% endif %}

{% if attributes %}
Property Details
----------------

{% for name in attributes %}
.. autoattribute:: {{ objname }}.{{ name }}

{% endfor %}
{% endif %}

{% if methods %}
Method Details
--------------

{% for name in methods %}
.. automethod:: {{ objname }}.{{ name }}

{% endfor %}
{% endif %}