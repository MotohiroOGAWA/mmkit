{{ fullname }}
{{ underline }}

{{ docstring }}

Properties
----------

.. autosummary::
   :nosignatures:

{% for name in attributes %}
   {{ fullname }}.{{ name }}
{% endfor %}

Methods
-------

.. autosummary::
   :nosignatures:

{% for name in methods %}
   {{ fullname }}.{{ name }}
{% endfor %}