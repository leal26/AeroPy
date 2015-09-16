.. pyFOIL Documentation

.. include xfoil_functions.rsf

Welcome to pyFOIL's documentation, an easy to use aerodynamic tool
==================================================================

Python interface for XFOIL. THe objective of this library is to be able to be
able to call XFOIL from Python iteratively for any simulation in a total of 4
lines total (one line for most uses). For a thorough explanation and tutorials
please check the table of contents.

Contents:

.. toctree::
   :numbered:
   :titlesonly:
   :maxdepth: 2

   xfoil_module
   aero_module
   main_module
   tutorial

To Do
======

- Include asymmetric wing
- Create airfoil and wing classes

Recommended Collaborations
==========================

Please use and adapt this library to your needs. There are several
functionalities I wished to implement, but did not have the time. Hence
I am *strongly* recommending the following collaborations:

- Airfoil generator with a GUI
- Atmospheric module (use the library already available in aero_module)
- Extend aero_module for wings with non-constant cross sections

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Examples
==================

.. code-block:: python

    from xfoil_tools import *
    find_coefficients(airfoil='naca0012', alpha=12.,NACA=True)
    >>{'CD',}
