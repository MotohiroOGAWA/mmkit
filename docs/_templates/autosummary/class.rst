{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :no-index-entry:

{% if docstring %}
{{ docstring }}
{% endif %}

{% if attributes %}
Properties
----------

.. autosummary::
   :nosignatures:

{% for name in attributes %}
   {{ fullname }}.{{ name }}
{% endfor %}
{% endif %}

{% if methods %}
Methods
-------

.. autosummary::
   :nosignatures:

{% for name in methods %}
   {{ fullname }}.{{ name }}
{% endfor %}
{% endif %}