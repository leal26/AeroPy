.. AeroPy Documentation

.. include xfoil_functions.rsf

Welcome to AeroPy's, an easy to use aerodynamic tool
==================================================================
AeroPy is an library for calculating aerodynamic properties. The main feature of
this library is the Python interface with XFOIL. The main objective of this library is to be able 
to use XFOIL via Python iteratively in a total of 4 lines total (one line for most uses). 
Through this interface coupling with other softwares (Abaqus, Ansys, etc) is possible
and iterative processes (optimization, design sensitivity) are possible.
For a thorough explanation please check the documentation and the tutorials.

Contents:

.. toctree::
   :numbered:
   :titlesonly:
   :maxdepth: 2

   main_module
   xfoil_module
   aero_module
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
