.. _dimtutorial:

=======================
Dimensionality analysis
=======================

This is a example of analysis of the dimensionality of a structure using
the :func:`ase.geometry.analyze_dimensionality` function. This is
useful for finding low-dimensional materials, such as 1D chain-like
structures, 2D layered structures, or structures with multiple dimensionality
types, such as 1D+3D.

The example below creates a layered :mol:`MoS_2` structure and analyzes its
dimensionality.

.. literalinclude:: dimexample.py

Coloring the atoms by their tags shows the distinct bonded clusters, which in
this case are separate layers.

The method is described in the article:

  | P.M. Larsen, M. Pandey, M. Strange, and K. W. Jacobsen
  | `Definition of a scoring parameter to identify low-dimensional materials components`__
  | 2018

__ http://arxiv.org/abs/1808.02114
