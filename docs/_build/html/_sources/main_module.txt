AeroPy Documentation
*********************

This project was elaborated because of the need of optimizing an airfoil
according to its aerodynamic and structural performance (conference paper).
Contrary to available options such as XFLR5, AeroPy aims to be an open source
Python code. It is a quick and easy way to find the aerodynamic pressure and the
drag and lift coefficients through the use of MIT's XFOIL embedded in Python.
Contrary to the other libraries, AeroPy does not have a GUI and is intended
to be used in Python codes for optimizations or any other process that requires
several aerodynamic analysis.

AeroPy is separated in three modules: (it is programmed in such a way that a new
module can be easily inserted)

- xfoil: contains all the functions relating to XFOIL.
- aero: contains all functions related to aerodynamics, but not related to XFOIL.
- PyXFOIL: imports all functions and acts as the intersection of all different
   modules

Code
====

.. automodule:: AeroPy
   :members:
