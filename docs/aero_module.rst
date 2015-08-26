Aerodynamic Module Documentation
***************************************

Lifting Line Theory
===================
If :math:`\theta_0` is an arbitrary span-wise location:

.. math::
   \alpha (\theta_o)=\frac{2b}{\pi c(\theta_o)} \sum_1^N A_n sin(n \theta_o) + \alpha_{L=0}(\theta_o) + \sum_1^N n A_n \frac{sin(n\theta_o)}{sin(\theta_o)}
   :label: LLT_full

Each equation has :math:`N` unknowns (:math:`A_n`), so if there are N :math:`\theta_o`, we have NxN system, which in Einstein notation can be written as:

.. math::
   C_{ij}A_{i}=D_{i}
   :label: LLT_simple

where, :math:`i=0,...,N`, :math:`j=0,...,N` and :

.. math::
   C_{ij}= \left( \frac{2b}{\pi c(j)} + \frac{n}{sin \theta(i)} \right) sin(n \theta(i))
   :label: C

.. math::
   A_i=A(i)
   :label: A

.. math::
   D_i=\alpha(i)-\alpha_{L=0}(i)
   :label: D

where  :math:`n=1,3,5,...,N-1`. Since we are considering a symmetric wing, all of the even terms would cancel each other

.. figure:: images/elliptical_LLT.png
   :align: center

The code
========
.. automodule:: aero_module
   :members:
